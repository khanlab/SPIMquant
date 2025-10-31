import pandas as pd

# read the original csv
df = pd.read_csv(snakemake.input.csv)

# identify rows with same left/right label
same_label = df[df["left label"] == df["right label"]]
diff_label = df[df["left label"] != df["right label"]]

# create left/right entries only for the differing-label rows
df_left = diff_label.copy()
df_left["index"] = df_left["left label"]
df_left["name"] = "left " + df_left["Structure"]
df_left["ABI"] = "left " + df_left["ABI"]

df_right = diff_label.copy()
df_right["index"] = df_right["right label"]
df_right["name"] = "right " + df_right["Structure"]
df_right["ABI"] = "right " + df_right["ABI"]

# rows where left/right labels are the same → use one entry
df_same = same_label.copy()
df_same["index"] = df_same["left label"]
df_same["name"] = df_same["Structure"]

# combine all
df_out = pd.concat([df_left, df_right, df_same], ignore_index=True)

# keep only the desired columns, rename for clarity
df_out = df_out[["index", "name", "hierarchy", "tissue type", "ABI"]]
df_out = df_out.rename(columns={"tissue type": "tissue_type"})

# sort by index
# df_out = df_out.sort_values(["index"]).reset_index(drop=True)

# save as TSV
df_out.to_csv(snakemake.output.tsv, sep="\t", index=False)
