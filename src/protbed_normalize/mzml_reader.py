"""Convert mzML/mzXML files into a peak-intensity table."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pyteomics import mzml, mzxml

from protbed_normalize.schema import Category, NormalizationMeta, NormalizedTable


def _centroid_peak_indices(mz: np.ndarray, intensity: np.ndarray) -> np.ndarray:
    """Return indices of local maxima above zero intensity."""
    if len(intensity) < 3:
        return np.where(intensity > 0)[0]
    greater_than_prev = intensity[1:-1] > intensity[:-2]
    greater_than_next = intensity[1:-1] > intensity[2:]
    positive = intensity[1:-1] > 0
    peaks = np.where(greater_than_prev & greater_than_next & positive)[0] + 1
    return peaks


def _bin_peaks(
    mz_values: list[np.ndarray],
    intensity_values: list[np.ndarray],
    ppm: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge peaks across spectra by m/z within a ppm tolerance.

    Returns unique m/z bins and summed intensities per spectrum.
    """
    # Collect all peaks.
    all_mz = np.concatenate(mz_values) if mz_values else np.array([])
    all_int = np.concatenate(intensity_values) if intensity_values else np.array([])
    if len(all_mz) == 0:
        return np.array([]), np.array([])

    order = np.argsort(all_mz)
    sorted_mz = all_mz[order]
    sorted_int = all_int[order]

    bin_mz = [sorted_mz[0]]
    bin_int = [sorted_int[0]]
    for mz, intensity in zip(sorted_mz[1:], sorted_int[1:]):
        tol = bin_mz[-1] * ppm / 1e6
        if abs(mz - bin_mz[-1]) <= tol:
            # Weighted average m/z, summed intensity.
            total = bin_int[-1] + intensity
            if total > 0:
                bin_mz[-1] = (bin_mz[-1] * bin_int[-1] + mz * intensity) / total
            bin_int[-1] = total
        else:
            bin_mz.append(mz)
            bin_int.append(intensity)
    return np.array(bin_mz), np.array(bin_int)


def read_mzml(
    path: str | Path,
    category: Category | str = Category.METABOLOMICS,
    aggregate: str = "sum",
    ppm: float = 10.0,
) -> NormalizedTable:
    """Read an mzML or mzXML file and produce a feature × sample table.

    Parameters
    ----------
    path: str | Path
        Path to .mzML or .mzXML file.
    category: Category | str
        Omics category label for metadata.
    aggregate: str
        How to aggregate intensities per spectrum: "sum" or "max".
    ppm: float
        Mass tolerance in parts per million for peak binning.

    Returns
    -------
    NormalizedTable
        A table where each row is a spectrum (sample) and columns are binned m/z peaks.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    category = Category(category)

    reader_cls = mzml.MzML if suffix == ".mzml" else mzxml.MzXML
    spectra: list[dict[str, Any]] = []
    with reader_cls(str(path)) as reader:
        for spec in reader:
            spectra.append(spec)

    if not spectra:
        raise ValueError(f"No spectra found in {path}")

    # Extract centroid peaks per spectrum.
    spectrum_peaks: list[tuple[str, np.ndarray, np.ndarray]] = []
    for spec in spectra:
        scan_id = spec.get("id", "spectrum")
        mz = np.array(spec["m/z array"], dtype=float)
        intensity = np.array(spec["intensity array"], dtype=float)
        peaks = _centroid_peak_indices(mz, intensity)
        spectrum_peaks.append((scan_id, mz[peaks], intensity[peaks]))

    # Build a common m/z grid from all peaks.
    all_mz = [peaks[1] for peaks in spectrum_peaks]
    all_int = [peaks[2] for peaks in spectrum_peaks]
    common_mz, _ = _bin_peaks(all_mz, all_int, ppm=ppm)

    if len(common_mz) == 0:
        raise ValueError("No peaks could be extracted from the file.")

    # Map each spectrum's peaks onto the common grid.
    rows = []
    for scan_id, mz, intensity in spectrum_peaks:
        row = np.zeros(len(common_mz), dtype=float)
        for m, i in zip(mz, intensity):
            # Find closest bin within tolerance.
            tol = m * ppm / 1e6
            diffs = np.abs(common_mz - m)
            idx = np.argmin(diffs)
            if diffs[idx] <= tol:
                if aggregate == "sum":
                    row[idx] += i
                elif aggregate == "max":
                    row[idx] = max(row[idx], i)
                else:
                    raise ValueError(f"Unknown aggregate: {aggregate}")
        rows.append((scan_id, row))

    df = pd.DataFrame(
        [r[1] for r in rows],
        columns=[f"mz_{m:.4f}" for m in common_mz],
    )
    df.insert(0, "spectrum_id", [r[0] for r in rows])

    meta = NormalizationMeta(
        source_path=str(path),
        detected_category=category,
        detected_orientation="mzml_spectra_as_rows",
        id_column="spectrum_id",
        sample_columns=df.columns.tolist()[1:],
        missing_strategy="none",
        duplicate_strategy=aggregate,
        warnings=[
            f"Extracted {len(spectra)} spectra and {len(common_mz)} unique m/z bins at {ppm} ppm.",
        ],
    )

    return NormalizedTable(data=df, meta=meta)
