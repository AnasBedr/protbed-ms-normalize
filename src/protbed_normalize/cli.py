"""Command-line interface for protbed-ms-normalize."""

from __future__ import annotations

import json
import os
from pathlib import Path

import click

from protbed_normalize.mzml_reader import read_mzml
from protbed_normalize.parser import normalize_table
from protbed_normalize.schema import Category


@click.command()
@click.argument("input_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file path. If omitted, written next to input with .normalized.csv suffix.",
)
@click.option(
    "-c",
    "--category",
    type=click.Choice([c.value for c in Category], case_sensitive=False),
    default=None,
    help="Omics category. Auto-detected if not provided.",
)
@click.option(
    "--duplicate-strategy",
    type=click.Choice(["mean", "median", "max", "first"], case_sensitive=False),
    default="mean",
    help="How to collapse duplicate feature IDs.",
)
@click.option(
    "--missing-strategy",
    type=click.Choice(["none", "drop", "zero", "mean", "median"], case_sensitive=False),
    default="none",
    help="How to handle missing intensity values.",
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["csv", "excel", "json"], case_sensitive=False),
    default=None,
    help="Output format. Inferred from output path if omitted.",
)
@click.option("--strip-prefix", is_flag=True, help="Remove the leading text shared by all sample names.")
@click.option("--alias-samples", is_flag=True, help="Rename samples to S1, S2, ... and save a mapping JSON.")
def main(
    input_path: Path,
    output: Path | None,
    category: str | None,
    duplicate_strategy: str,
    missing_strategy: str,
    fmt: str | None,
    strip_prefix: bool,
    alias_samples: bool,
) -> None:
    """Normalize a mass spectrometry table or mzML file for downstream analysis."""
    suffix = input_path.suffix.lower()
    cat = Category(category) if category else None

    if suffix in {".mzml", ".mzxml"}:
        table = read_mzml(input_path, category=cat or Category.METABOLOMICS)
    else:
        table = normalize_table(
            input_path,
            category=cat,
            duplicate_strategy=duplicate_strategy,
            missing_strategy=missing_strategy,
        )

    if output is None:
        output = input_path.with_suffix(".normalized.csv")

    mapping: dict[str, str] = {}
    samples = list(table.meta.sample_columns)
    if (strip_prefix or alias_samples) and samples:
        if alias_samples:
            new = [f"S{i + 1}" for i in range(len(samples))]
        else:
            prefix = os.path.commonprefix(samples) if len(samples) > 1 else ""
            cut = max(prefix.rfind("_"), prefix.rfind("-"), prefix.rfind(" ")) + 1
            new = [s[cut:] or s for s in samples]
            if len(set(new)) != len(new):
                new = samples
        mapping = {n: o for n, o in zip(new, samples) if n != o}
        display = {n: (o[4:] if o.startswith("col_") and o[4:5].isdigit() else o) for n, o in mapping.items()}
        table.data = table.data.rename(columns={o: n for n, o in mapping.items()})
        table.meta.sample_columns = new

    table.write(output, format=fmt)
    if alias_samples and mapping:
        map_path = output.with_suffix(".sample_mapping.json")
        map_path.write_text(json.dumps(display, indent=2), encoding="utf-8")
        click.echo(f"Sample mapping written to: {map_path}")

    click.echo(f"Normalized table written to: {output}")
    click.echo(f"Detected category: {table.meta.detected_category.value}")
    click.echo(f"Detected orientation: {table.meta.detected_orientation}")
    click.echo(f"Identifier column: {table.meta.id_column}")
    click.echo(f"Sample columns: {len(table.meta.sample_columns)}")
    if table.meta.warnings:
        click.echo("Warnings:")
        for warning in table.meta.warnings:
            click.echo(f"  - {warning}")
    click.echo(
        "\nFor downstream statistics, visualization, and AI interpretation, "
        "upload the normalized file to ProtBed: https://protbed.com"
    )


if __name__ == "__main__":
    main()
