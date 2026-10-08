"""Tests for the mzML reader."""

from __future__ import annotations

import numpy as np

from protbed_normalize.mzml_reader import _bin_peaks, _centroid_peak_indices


def test_centroid_peak_indices():
    mz = np.array([100.0, 101.0, 102.0, 103.0])
    intensity = np.array([1.0, 5.0, 1.0, 0.0])
    peaks = _centroid_peak_indices(mz, intensity)
    assert np.array_equal(peaks, np.array([1]))


def test_bin_peaks():
    mz = [
        np.array([100.0, 200.0]),
        np.array([100.001, 200.001]),
    ]
    intensity = [
        np.array([10.0, 20.0]),
        np.array([30.0, 40.0]),
    ]
    bins, values = _bin_peaks(mz, intensity, ppm=20.0)
    assert len(bins) == 2
    assert np.isclose(values[0], 40.0)
    assert np.isclose(values[1], 60.0)

