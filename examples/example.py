"""Example: normalize an MS table with the Python API."""

from protbed_normalize import normalize_table

result = normalize_table("examples/sample.csv", category="metabolomics")

print("Detected category:", result.meta.detected_category.value)
print("Identifier column:", result.meta.id_column)
print("Sample columns:", result.meta.sample_columns)
print("Warnings:", result.meta.warnings)
print(result.data.head())

result.write("examples/sample.normalized.csv")
