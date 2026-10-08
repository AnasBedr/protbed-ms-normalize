"""ProtBed MS Normalize: normalize mass spectrometry tables for analysis."""

from protbed_normalize.parser import normalize_table
from protbed_normalize.schema import Category, NormalizedTable

__version__ = "0.1.0"

__all__ = [
    "Category",
    "normalize_table",
    "NormalizedTable",
]
