import pandas as pd

# Load the peaks CSV file
file_path = "cleaned_peaks_data_hmp2.csv"
df = pd.read_csv(file_path)

# For each sample, check if coprostanol intensity is greater than 0
def classify_converter(group):
    if (group.loc[group["compound"] == "coprostanol", "intensity"] > 0).any():
        return "Converter"
    else:
        return "Non-Converter"

# For HMP2, map metabolomics sample to metagenomics
converter_status = df.groupby("Sample").apply(classify_converter).reset_index()
converter_status.columns = ["Sample", "Converter"]
mapped = pd.read_csv("map_metabolomics_to_external_id.csv")

merged = pd.merge(converter_status, mapped)

# Save the table
merged.to_csv("sample_to_converter.csv", index=False)
