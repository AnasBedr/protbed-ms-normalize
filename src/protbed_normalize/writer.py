"""Write normalized tables to common formats."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from protbed_normalize.schema import NormalizedTable


def _infer_fmt(path: str | Path) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return "csv"
    if suffix in {".xlsx", ".xls"}:
        return "excel"
    if suffix == ".json":
        return "json"
    raise ValueError(f"Cannot infer output format from path: {path}")


def write_table(
    table: "NormalizedTable",
    output_path: str | Path,
    fmt: str | None = None,
) -> None:
    """Write a NormalizedTable to disk.

    CSV and Excel write the data matrix. JSON writes data as records plus metadata.
    A sidecar `.meta.json` is always written for CSV/Excel outputs.
    """
    path = Path(output_path)
    fmt = (fmt or _infer_fmt(path)).lower()

    if fmt == "csv":
        table.data.to_csv(path, index=False)
        _write_meta(table, path.with_suffix(".meta.json"))
    elif fmt == "excel":
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            table.data.to_excel(writer, sheet_name="normalized", index=False)
        _write_meta(table, path.with_suffix(".meta.json"))
    elif fmt == "json":
        payload = {
            "meta": table.meta.model_dump(),
            "data": table.data.to_dict(orient="records"),
        }
        path.write_text(json.dumps(payload, indent=2))
    else:
        raise ValueError(f"Unsupported output format: {fmt}")


def _write_meta(table: "NormalizedTable", path: Path) -> None:
    path.write_text(json.dumps(table.meta.model_dump(), indent=2))
