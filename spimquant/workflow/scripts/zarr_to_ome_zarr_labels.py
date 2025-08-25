import json
import zarr
import dask.array as da
import pandas as pd
from ome_zarr.writer import write_labels, write_label_metadata
from ome_zarr.scale import Scaler
from dask.diagnostics import ProgressBar


in_zarr = snakemake.input.zarr
metadata_json = snakemake.input.metadata_json
label_tsv = snakemake.input.label_tsv
max_layer = (
    snakemake.params.max_downsampling_layers
)  # number of downsamplings by 2 to include in zarr
rechunk_size = snakemake.params.rechunk_size
scaling_method = snakemake.params.scaling_method
label_name = snakemake.params.label_name
out_zarr = snakemake.output.zarr


darr = da.from_zarr(in_zarr)

# prepare metadata for ome-zarr
with open(metadata_json) as fp:
    metadata = json.load(fp)

voxdim = [
    float(metadata["physical_size_z"]) / 1000.0,
    float(metadata["physical_size_y"]) / 1000.0,
    float(metadata["physical_size_x"]) / 1000.0,
]


coordinate_transformations = []
# for each resolution (dataset), we have a list of dicts, transformations to apply..
# in this case just a single one (scaling by voxel size)

for l in range(max_layer + 1):

    coordinate_transformations.append(
        [
            {
                "scale": [
                    voxdim[0],
                    (2**l) * voxdim[1],
                    (2**l) * voxdim[2],
                ],  # image-pyramids in XY only
                "type": "scale",
            }
        ]
    )


axes = [{"name": ax, "type": "space", "unit": "micrometer"} for ax in ["z", "y", "x"]]


store = zarr.DirectoryStore(out_zarr)
root = zarr.group(store, path="/", overwrite=True)
scaler = Scaler(max_layer=max_layer, method="nearest")


# label metadata - convert from bids tsv to list of dicts
df = pd.read_csv(label_tsv, sep="\t", index_col="index")

colors = []
properties = []
for i, row in df.iterrows():
    if i > 0:  # leave out bg label
        hex = row["color"].lstrip("#")
        rgba = tuple(int(hex[i : i + 2], 16) for i in (0, 2, 4))
        colors.append({"label-value": i, "rgba": rgba})
        properties.append(
            {"label-value": i, "name": row["name"], "abbreviation": row["abbreviation"]}
        )


with ProgressBar():
    write_labels(
        labels=darr,
        group=root,
        scaler=scaler,
        name=label_name,
        coordinate_transformations=coordinate_transformations,
        axes=axes,
    )

    write_label_metadata(
        group=root, name=f"/labels/{label_name}", colors=colors, properties=properties
    )
