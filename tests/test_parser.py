"""Tests for the tabular parser."""

from __future__ import annotations

import pandas as pd

from protbed_normalize.parser import normalize_table
from protbed_normalize.schema import Category


def test_normalize_compounds_as_rows(tmp_path):
    df = pd.DataFrame({
        "Compound Name": ["A", "B", "C"],
        "m/z": ["100.1", "200.2", "300.3"],
        "RT": ["1.1", "2.2", "3.3"],
        "Sample_1": ["10", "20", "30"],
        "Sample_2": ["15", "25", "35"],
    })
    path = tmp_path / "input.csv"
    df.to_csv(path, index=False)

    result = normalize_table(path, category=Category.METABOLOMICS)

    assert result.meta.detected_category == Category.METABOLOMICS
    assert result.meta.detected_orientation == "compounds_as_rows"
    assert result.meta.id_column in result.data.columns
    assert len(result.meta.sample_columns) == 2
    assert result.data[result.meta.sample_columns[0]].dtype.kind == "f"


def test_normalize_samples_as_rows(tmp_path):
    df = pd.DataFrame({
        "sample": ["Sample_1", "Sample_2", "Sample_3"],
        "A": ["10", "15", "12"],
        "B": ["20", "25", "22"],
        "C": ["30", "35", "32"],
    })
    path = tmp_path / "input.csv"
    df.to_csv(path, index=False)

    result = normalize_table(path)

    assert result.meta.detected_orientation == "samples_as_rows"
    # After transposition, original row labels become sample columns.
    assert any("Sample_1" in c for c in result.data.columns)
    # Original feature labels move into the identifier column values.
    id_values = result.data[result.meta.id_column].astype(str).tolist()
    assert "A" in id_values
    assert "B" in id_values
    assert "C" in id_values



def test_duplicate_collapse(tmp_path):
    df = pd.DataFrame({
        "Name": ["A", "A", "B"],
        "Sample_1": ["10", "20", "30"],
        "Sample_2": ["5", "15", "25"],
    })
    path = tmp_path / "input.csv"
    df.to_csv(path, index=False)

    result = normalize_table(path, duplicate_strategy="mean")

    assert result.meta.duplicate_count == 1
    row = result.data[result.data[result.meta.id_column] == "A"]
    assert row["Sample_1"].iloc[0] == 15.0
    assert row["Sample_2"].iloc[0] == 10.0


def test_missing_imputation(tmp_path):
    df = pd.DataFrame({
        "Name": ["A", "B", "C"],
        "Sample_1": ["10", "", "30"],
        "Sample_2": ["5", "15", ""],
    })
    path = tmp_path / "input.csv"
    df.to_csv(path, index=False)

    result = normalize_table(path, missing_strategy="mean")

    assert not result.data[result.meta.sample_columns].isna().any().any()
