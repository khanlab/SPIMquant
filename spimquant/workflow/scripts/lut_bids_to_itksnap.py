import json

import pandas as pd

lut_header = """
################################################
# ITK-SnAP Label Description File
# File format:
# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL
# Fields:
#    IDX:   Zero-based index
#    -R-:   Red color component (0..255)
#    -G-:   Green color component (0..255)
#    -B-:   Blue color component (0..255)
#    -A-:   Label transparency (0.00 .. 1.00)
#    VIS:   Label visibility (0 or 1)
#    IDX:   Label mesh visibility (0 or 1)
#  LABEL:   Label description
################################################
"""


# Function to convert hex RGB string to 8-bit intensity values
def hex_to_rgb(hex_string):
    r = int(hex_string[0:2], 16)
    g = int(hex_string[2:4], 16)
    b = int(hex_string[4:6], 16)
    return r, g, b


df = pd.read_csv(snakemake.input.tsv, sep="\t")

# Apply the function to create new 'R', 'G', and 'B' columns
df[["R", "G", "B"]] = df["color"].apply(lambda x: pd.Series(hex_to_rgb(x)))

df["A"] = 1
df["VIS"] = 1
df["IDX"] = 1
df["LABEL"] = df["name"].apply(lambda x: f'"{x}"')
df.iloc[0]["A"] = 0
df.iloc[0]["VIS"] = 0
df.iloc[0]["IDX"] = 0


# write header first
with open(snakemake.output.lut, "w") as f:
    f.write(lut_header)


df[["index", "R", "G", "B", "A", "VIS", "IDX", "LABEL"]].to_csv(
    snakemake.output.lut, sep="\t", index=False, header=False, mode="a", quoting=3
)
