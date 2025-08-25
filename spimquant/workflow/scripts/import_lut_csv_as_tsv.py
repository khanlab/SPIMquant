import pandas as pd

df = pd.read_csv(
    snakemake.input.csv, header=None, names=["abbreviation", "name", "index"]
)
df.to_csv(
    snakemake.output.tsv,
    sep="\t",
    index=False,
    columns=["index", "name", "abbreviation"],
)
