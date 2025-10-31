import pandas as pd
import random
from matplotlib import colormaps
from matplotlib.colors import to_hex

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


def hex_to_rgb(hex_string):
    try:
        hex_string = str(hex_string).strip().replace("#", "")
        if len(hex_string) != 6:
            raise ValueError
        r = int(hex_string[0:2], 16)
        g = int(hex_string[2:4], 16)
        b = int(hex_string[4:6], 16)
        return r, g, b
    except Exception:
        return None


# Read TSV
df = pd.read_csv(snakemake.input.tsv, sep="\t")

# Build fallback colors using tab20 colormap, sampled evenly
colormap = colormaps.get_cmap("tab20")
n = len(df)
fallback_colors = [to_hex(colormap(i / max(1, n - 1))) for i in range(n)]

valid_colors = []
for i, row in df.iterrows():
    color = row.get("color", "")
    rgb = hex_to_rgb(color)
    if rgb is None:
        # fallback: use colormap or random color if needed
        fallback_hex = (
            fallback_colors[i]
            if i < len(fallback_colors)
            else "#{:06x}".format(random.randint(0, 0xFFFFFF))
        )
        rgb = hex_to_rgb(fallback_hex.replace("#", ""))
    valid_colors.append(rgb)

# Unpack into separate columns
df["R"], df["G"], df["B"] = zip(*valid_colors)

# Add other columns
df["A"] = 1
df["VIS"] = 1
df["IDX"] = 1
df["LABEL"] = df["name"].apply(lambda x: f'"{x}"')

# Make background (first row) invisible
df.iloc[0, df.columns.get_loc("A")] = 0
df.iloc[0, df.columns.get_loc("VIS")] = 0
df.iloc[0, df.columns.get_loc("IDX")] = 0

# Write output
with open(snakemake.output.lut, "w") as f:
    f.write(lut_header)

df[["index", "R", "G", "B", "A", "VIS", "IDX", "LABEL"]].to_csv(
    snakemake.output.lut, sep="\t", index=False, header=False, mode="a", quoting=3
)
