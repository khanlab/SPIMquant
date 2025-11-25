# import ngff_zarr as nz
from zarrnii import ZarrNii

# Read from directory store
# multiscales = nz.from_ngff_zarr(snakemake.input.zarr)

# print(multiscales)
# Write as .ozx file
# nz.to_ngff_zarr(snakemake.output.ozx, multiscales, version='0.5')

# note, this recomputes the multiscales
ZarrNii.from_ome_zarr(snakemake.input.zarr).to_ome_zarr(snakemake.output.ozx)
