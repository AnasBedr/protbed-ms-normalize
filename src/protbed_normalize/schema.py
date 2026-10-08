"""Schema, enums, and detection keywords for normalized MS tables."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class Category(str, Enum):
    """Supported omics categories."""

    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    LIPIDOMICS = "lipidomics"
    EXPOSOMICS = "exposomics"
    UNKNOWN = "unknown"


# Detection keywords are lower-cased before matching.
ID_KEYWORDS: dict[Category, list[list[str]]] = {
    Category.PROTEOMICS: [
        ["gene symbol", "gene name", "gene"],
        ["accession", "uniprot", "protein id", "entry", "protein"],
    ],
    Category.METABOLOMICS: [
        ["name", "compound", "metabolite", "feature"],
    ],
    Category.LIPIDOMICS: [
        ["name", "lipid", "species", "feature"],
    ],
    Category.EXPOSOMICS: [
        ["name", "compound", "chemical", "feature"],
    ],
    Category.UNKNOWN: [
        ["name", "id", "feature"],
    ],
}

MZ_KEYWORDS = [
    "m/z",
    "mz",
    "mass",
    "mass-to-charge",
    "precursor mz",
    "precursor_mz",
]

RT_KEYWORDS = [
    "rt",
    "retention time",
    "retention_time",
    "retention",
    "rtime",
]

SAMPLE_BLOCKLIST = [
    "blank",
    "qc",
    "pool",
    "standard",
    "std",
]


class NormalizationMeta(BaseModel):
    """Metadata describing how a table was normalized."""

    source_path: str
    detected_category: Category = Category.UNKNOWN
    detected_orientation: str = "unknown"
    id_column: str | None = None
    mz_column: str | None = None
    rt_column: str | None = None
    sample_columns: list[str] = Field(default_factory=list)
    duplicate_count: int = 0
    duplicate_strategy: str = "mean"
    missing_strategy: str = "none"
    warnings: list[str] = Field(default_factory=list)


class NormalizedTable(BaseModel):
    """A normalized MS table ready for downstream analysis."""

    data: Any  # pandas DataFrame
    meta: NormalizationMeta

    model_config = {"arbitrary_types_allowed": True}


    def write(
        self,
        output_path: str | Path,
        format: str | None = None,
    ) -> None:
        """Write the normalized table and its metadata to disk.

        Parameters
        ----------
        output_path: str | Path
            Destination path. Supported extensions: .csv, .xlsx, .json.
        format: str | None
            Optional override ("csv", "excel", "json"). If None, inferred from path.
        """
        from protbed_normalize.writer import write_table

        write_table(self, output_path, fmt=format)
