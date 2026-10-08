"""Parse and normalize messy MS tabular exports."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from protbed_normalize.schema import (
    ID_KEYWORDS,
    MZ_KEYWORDS,
    RT_KEYWORDS,
    SAMPLE_BLOCKLIST,
    Category,
    NormalizationMeta,
    NormalizedTable,
)


def _lower_strip(name: str) -> str:
    return name.lower().strip().replace("_", " ")


def _matches_keyword(name: str, keyword: str) -> bool:
    """Check whether a column name matches a keyword phrase."""
    normalized = _lower_strip(name)
    keyword_norm = _lower_strip(keyword)
    # Exact match or word-boundary match.
    if normalized == keyword_norm:
        return True
    if re.search(r"\b" + re.escape(keyword_norm) + r"\b", normalized):
        return True
    return False


def _find_column(df: pd.DataFrame, keywords: list[str]) -> str | None:
    """Return the first column name that matches any keyword."""
    for col in df.columns:
        for kw in keywords:
            if _matches_keyword(col, kw):
                return col
    return None


def _find_id_column(df: pd.DataFrame, category: Category) -> str | None:
    """Find the best identifier column for the given category."""
    keyword_groups = ID_KEYWORDS.get(category, ID_KEYWORDS[Category.UNKNOWN])
    for group in keyword_groups:
        match = _find_column(df, group)
        if match is not None:
            return match
    return None


def _detect_orientation(df: pd.DataFrame) -> str:
    """Heuristic: are rows compounds (features) or samples?"""
    n_rows, n_cols = df.shape
    sample_keywords = ["sample", "rep", "qc", "blank", "control"]

    # Check the index (row labels) for sample-like identifiers.
    row_id_score = sum(
        1
        for idx in df.index[: min(20, n_rows)]
        if any(kw in str(idx).lower() for kw in sample_keywords)
    )
    # Check column headers for sample-like identifiers.
    col_id_score = sum(
        1
        for col in df.columns[: min(20, n_cols)]
        if any(kw in str(col).lower() for kw in sample_keywords)
    )
    # Check the first column's values; if they look like sample IDs, rows are likely samples.
    first_col = df.columns[0]
    first_col_score = sum(
        1
        for val in df[first_col].dropna().head(20)
        if any(kw in str(val).lower() for kw in sample_keywords)
    )

    if (row_id_score + first_col_score) > col_id_score:
        return "samples_as_rows"
    return "compounds_as_rows"



def _looks_like_sample_column(
    name: str,
    id_col: str | None,
    mz_col: str | None,
    rt_col: str | None,
) -> bool:
    """Return True if a column is likely a sample intensity column."""
    lower = name.lower()
    # Exclude known metadata/ID columns.
    if id_col and name == id_col:
        return False
    if mz_col and name == mz_col:
        return False
    if rt_col and name == rt_col:
        return False
    # Exclude common metadata columns.
    metadata_keywords = [
        "mz", "rt", "retention", "name", "id", "accession",
        "gene", "protein", "compound", "formula", "adduct",
    ]
    tokens = set(re.split(r"[^0-9a-z]+", lower))
    if tokens & set(metadata_keywords) or "m/z" in lower:
        return False

    # Exclude QC/blank columns from default sample set.
    if any(kw in lower for kw in SAMPLE_BLOCKLIST):
        return False
    return True


def _reconcile_duplicates(
    df: pd.DataFrame,
    id_col: str,
    strategy: str,
) -> tuple[pd.DataFrame, int]:
    """Collapse duplicate feature IDs using the chosen strategy."""
    before = len(df)
    if id_col not in df.columns:
        return df, 0
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    non_numeric_cols = [c for c in df.columns if c not in numeric_cols]
    if id_col in numeric_cols:
        numeric_cols.remove(id_col)

    if strategy == "first":
        agg: dict[str, Any] = {c: "first" for c in non_numeric_cols}
        agg.update({c: "first" for c in numeric_cols})
    elif strategy == "mean":
        agg = {c: "first" for c in non_numeric_cols}
        agg.update({c: "mean" for c in numeric_cols})
    elif strategy == "median":
        agg = {c: "first" for c in non_numeric_cols}
        agg.update({c: "median" for c in numeric_cols})
    elif strategy == "max":
        agg = {c: "first" for c in non_numeric_cols}
        agg.update({c: "max" for c in numeric_cols})
    else:
        raise ValueError(f"Unknown duplicate strategy: {strategy}")

    grouped = df.groupby(id_col, as_index=False).agg(agg)
    after = len(grouped)
    return grouped, before - after


def _handle_missing_values(df: pd.DataFrame, strategy: str) -> pd.DataFrame:
    """Impute or drop missing values in numeric sample columns."""
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if strategy == "drop":
        df = df.dropna(subset=numeric_cols)
    elif strategy == "zero":
        df[numeric_cols] = df[numeric_cols].fillna(0)
    elif strategy == "mean":
        df[numeric_cols] = df[numeric_cols].fillna(df[numeric_cols].mean())
    elif strategy == "median":
        df[numeric_cols] = df[numeric_cols].fillna(df[numeric_cols].median())
    elif strategy == "none":
        pass
    else:
        raise ValueError(f"Unknown missing value strategy: {strategy}")
    return df


def _normalize_column_names(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    """Rename columns to shell-safe, canonical names where possible.

    Returns the renamed DataFrame and a mapping of original column names to new names.
    """
    rename_map: dict[str, str] = {}
    for col in df.columns:
        safe = re.sub(r"[^0-9a-zA-Z_]+", "_", str(col)).strip("_")
        if not safe:
            safe = "unnamed"
        if not safe[0].isalpha() and safe[0] != "_":
            safe = "col_" + safe
        # Avoid collisions.
        counter = 1
        base = safe
        existing = set(df.columns) | set(rename_map.values())
        while safe in existing and safe != col:
            safe = f"{base}_{counter}"
            counter += 1
        if safe != col:
            rename_map[col] = safe
    return df.rename(columns=rename_map), rename_map



def _read_input(path: str | Path) -> pd.DataFrame:
    """Read CSV or Excel into a DataFrame."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str)
    if suffix in {".xlsx", ".xls"}:
        # Read the first non-empty sheet.
        xl = pd.ExcelFile(path)
        for sheet in xl.sheet_names:
            df = pd.read_excel(xl, sheet_name=sheet, dtype=str)
            if not df.empty:
                return df
        raise ValueError("All Excel sheets are empty.")
    raise ValueError(f"Unsupported tabular input format: {suffix}")


def normalize_table(
    path: str | Path,
    category: Category | str | None = None,
    duplicate_strategy: str = "mean",
    missing_strategy: str = "none",
) -> NormalizedTable:
    """Normalize a tabular MS export into a clean compounds × samples matrix.

    Parameters
    ----------
    path: str | Path
        Path to a .csv or .xlsx file.
    category: Category | str | None
        Omics category. If None, inferred from column names.
    duplicate_strategy: str
        How to collapse duplicate feature IDs: "mean", "median", "max", or "first".
    missing_strategy: str
        How to handle missing values: "none", "drop", "zero", "mean", or "median".

    Returns
    -------
    NormalizedTable
        The normalized matrix plus metadata describing the transformations applied.
    """
    path = Path(path)
    df = _read_input(path)
    warnings: list[str] = []

    # Drop completely empty rows/columns.
    df = df.dropna(how="all").dropna(axis=1, how="all")

    # Infer category if not provided.
    if category is None:
        header_text = " ".join(str(c).lower() for c in df.columns)
        if any(kw in header_text for kw in ["accession", "uniprot", "gene symbol", "protein"]):
            category = Category.PROTEOMICS
        elif any(kw in header_text for kw in ["lipid", "lipidomics"]):
            category = Category.LIPIDOMICS
        elif any(kw in header_text for kw in ["exposome", "exposomics", "chemical"]):
            category = Category.EXPOSOMICS
        elif any(kw in header_text for kw in ["metabolite", "compound", "mz", "rt"]):
            category = Category.METABOLOMICS
        else:
            category = Category.UNKNOWN
            warnings.append("Could not infer omics category; using 'unknown'.")
    else:
        category = Category(category)

    # Detect and fix orientation.
    orientation = _detect_orientation(df)
    if orientation == "samples_as_rows":
        df = df.T.reset_index()
        df.columns = [str(c) for c in df.iloc[0]]
        df = df.iloc[1:].reset_index(drop=True)
        warnings.append("Detected samples-as-rows orientation and transposed the table.")

    # Find identifier and metadata columns.
    id_col = _find_id_column(df, category)
    mz_col = _find_column(df, MZ_KEYWORDS)
    rt_col = _find_column(df, RT_KEYWORDS)

    if id_col is None:
        # Fall back to the first column if nothing else matches.
        id_col = df.columns[0]
        warnings.append(f"No clear identifier column found; using first column '{id_col}'.")

    # Identify sample intensity columns.
    sample_cols = [
        c for c in df.columns
        if _looks_like_sample_column(c, id_col, mz_col, rt_col)
    ]

    if not sample_cols:
        warnings.append("No sample intensity columns were detected.")

    # Convert sample columns to numeric.
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Reconcile duplicates.
    df, dup_count = _reconcile_duplicates(df, id_col, duplicate_strategy)
    if dup_count:
        warnings.append(
            f"Collapsed {dup_count} duplicate feature IDs using '{duplicate_strategy}'."
        )


    # Handle missing values.
    df = _handle_missing_values(df, missing_strategy)
    if missing_strategy != "none":
        warnings.append(f"Applied '{missing_strategy}' imputation to missing values.")

    # Rename columns for safety and update metadata to match.
    df, rename_map = _normalize_column_names(df)
    id_col = rename_map.get(id_col, id_col) if id_col else None
    mz_col = rename_map.get(mz_col, mz_col) if mz_col else None
    rt_col = rename_map.get(rt_col, rt_col) if rt_col else None
    sample_cols = [rename_map.get(c, c) for c in sample_cols]


    meta = NormalizationMeta(
        source_path=str(path),
        detected_category=category,
        detected_orientation=orientation,
        id_column=id_col,
        mz_column=mz_col,
        rt_column=rt_col,
        sample_columns=sample_cols,
        duplicate_count=dup_count,
        duplicate_strategy=duplicate_strategy,
        missing_strategy=missing_strategy,
        warnings=warnings,
    )

    return NormalizedTable(data=df, meta=meta)
