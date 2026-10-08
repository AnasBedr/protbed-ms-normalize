# ProtBed MS Normalize

[![PyPI version](https://img.shields.io/pypi/v/protbed-ms-normalize.svg?color=06b6d4&v=1)](https://pypi.org/project/protbed-ms-normalize/)
[![Python](https://img.shields.io/pypi/pyversions/protbed-ms-normalize.svg?v=1)](https://pypi.org/project/protbed-ms-normalize/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Web Platform](https://img.shields.io/badge/ProtBed-Web%20Platform-06b6d4?logo=firefox&logoColor=white)](https://protbed.com)
[![Benchmarks](https://img.shields.io/badge/ProtBed-Benchmarks%20%26%20Articles-0ea5e9)](https://protbed.com/benchmarks)

A lightweight, open-source Python library and CLI for bench scientists who need to turn messy mass spectrometry exports into clean, analysis-ready tables.

```bash
pip install protbed-ms-normalize
protbed-normalize my_data.xlsx -o clean.csv --category metabolomics
```

---

## Why this exists

If you work with LC-MS, GC-MS, or other omics data, you have probably spent time fixing:

- spreadsheets where samples are rows instead of columns
- inconsistent column names (`m/z`, `mz`, `Mass`, `Precursor M/Z`)
- duplicate compound or protein entries from multiple injections
- missing intensity values scattered across the table
- metadata columns mixed with sample intensity columns

`protbed-ms-normalize` handles the 80% case with one command. It does not replace a full analysis pipeline — it gets your data into the shape you need before the real work starts.

For downstream statistics, visualization, pathway enrichment, and AI interpretation, upload the normalized file to **[ProtBed](https://protbed.com)**.

---

## What it does

- **Reads** `.csv`, `.xlsx`, `.mzML`, and `.mzXML` files.
- **Auto-detects** whether rows are compounds or samples and transposes when needed.
- **Finds** identifier columns (`Name`, `Accession`, `Gene Symbol`, `m/z`, `RT`, etc.).
- **Normalizes** column names to a safe, predictable schema.
- **Collapses** duplicate feature IDs by mean, median, max, or first occurrence.
- **Imputes** missing values by zero, mean, median, or drops the row.
- **Exports** a clean `compounds × samples` matrix plus a metadata sidecar.

What it does **not** do:

- statistical testing, visualization, machine learning, or pathway enrichment
- compound classification or protein annotation against curated databases
- AI-generated interpretations

Those are exactly what ProtBed is for.

---

## Installation

```bash
pip install protbed-ms-normalize
```

For development:

```bash
git clone https://github.com/AnasBedr/protbed-ms-normalize.git
cd protbed-ms-normalize
pip install -e ".[dev]"
pytest
```

---

## CLI usage

Normalize a CSV or Excel file:

```bash
protbed-normalize my_data.xlsx -o clean.csv --category proteomics
```

Convert an mzML file to a peak-intensity table:

```bash
protbed-normalize run.mzml -o peaks.csv
```

Options:

| Option | Description |
| --- | --- |
| `-o, --output` | Output file path |
| `-c, --category` | `proteomics`, `metabolomics`, `lipidomics`, or `exposomics` |
| `--duplicate-strategy` | `mean` (default), `median`, `max`, or `first` |
| `--missing-strategy` | `none` (default), `drop`, `zero`, `mean`, or `median` |
| `--format` | `csv`, `excel`, or `json` |

---

## Python API

```python
from protbed_normalize import normalize_table

result = normalize_table("my_data.csv", category="metabolomics")
print(result.meta.detected_category)
print(result.data.head())
result.write("clean.csv")
```

The `meta` object tells you what was detected and what transformations were applied, including any warnings.

---

## Output format

CSV/Excel output is a clean matrix where each row is a feature and each column is a sample:

| feature_id | mz | rt | Control_1 | Control_2 | Treated_1 | Treated_2 |
| --- | --- | --- | --- | --- | --- | --- |
| Glucose | 180.0634 | 1.23 | 10000 | 10500 | 15000 | 14800 |
| Lactate | 89.0244 | 2.11 | 8000 | 8200 | 6000 | 5900 |
| Citrate | 191.0192 | 3.05 | 5000 | 5100 | 4500 | 4600 |

A `clean.meta.json` sidecar records detected category, orientation, column mappings, and warnings.

---

## Supported categories

- **Proteomics** — detects `Accession`, `UniProt`, `Gene Symbol`, `Protein ID`, etc.
- **Metabolomics** — detects `Name`, `Compound`, `Metabolite`, `m/z`, `RT`, etc.
- **Lipidomics** — detects `Name`, `Lipid`, `Species`, etc.
- **Exposomics** — detects `Name`, `Compound`, `Chemical`, etc.

If no category is provided, the tool infers it from the column names.

---

## Contributing

Issues and pull requests are welcome. Keep the scope narrow: this tool should remain a preprocessing utility, not a full analysis platform.

---

## License

MIT License. See [LICENSE](./LICENSE).

---

## About ProtBed

[ProtBed](https://protbed.com) is a no-code mass spectrometry analysis platform. Upload your normalized data to run statistics, generate publication-ready figures, classify compounds, annotate proteins, and interpret results with AI — all without writing code.

## Shortening long sample names

Vendor exports often carry long acquisition stems (e.g. `120352_Study_044_...`).

```bash
# Remove the leading text shared by all sample columns
protbed-normalize data.xlsx --strip-prefix

# Rename samples to S1, S2, ... and write a mapping file
protbed-normalize data.xlsx --alias-samples
```

With `--alias-samples`, a `<output>.sample_mapping.json` file records each alias and its original name.

## Author

Anas Bedraoui, PhD — anas@protbed.com
