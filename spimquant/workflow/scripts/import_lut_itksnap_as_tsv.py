#!/usr/bin/env python3
import pandas as pd


# Function to convert 8‐bit R, G, B values to a hex string (e.g., "ff0000")
def rgb_to_hex(r, g, b):
    return "{:02x}{:02x}{:02x}".format(int(r), int(g), int(b))


# Read the ITK‐Snap LUT file.
# The ITK‐Snap LUT contains a header (all lines starting with "#")
# so we use the "comment" argument to skip those lines.
# The data columns are (in order):
#   index, R, G, B, A, VIS, MSH, LABEL
columns = ["index", "R", "G", "B", "A", "VIS", "MSH", "LABEL"]
df = pd.read_csv(
    snakemake.input.lut, sep=r"\s+", comment="#", header=None, names=columns
)

# Remove the surrounding quotes from the LABEL column to get the ROI name.
df["name"] = df["LABEL"].str.strip('"')

# Create the BIDS color column by converting the R, G, B values to a hex string.
df["color"] = df.apply(lambda row: rgb_to_hex(row["R"], row["G"], row["B"]), axis=1)

# For the BIDS LUT, we only need the columns "index", "name", and "color"
df_bids = df[["index", "name", "color"]]

# Write out the BIDS dseg LUT file as a TSV.
df_bids.to_csv(snakemake.output.tsv, sep="\t", index=False)
