import json
import pandas as pd

col_order = ["name", "abbreviation", "color"]

with open(snakemake.input.json) as fp:
    lut = json.load(fp)

(
    pd.DataFrame(lut, columns=["index", "color", "abbreviation", "name"])
    .set_index("index")[col_order]
    .to_csv(snakemake.output.tsv, sep="\t")
)
