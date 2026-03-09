from ngff_zarr.rfc9_zip import write_store_to_zip
from zarr.storage import LocalStore

# Direct conversion of existing store to .ozx
source_store = LocalStore(snakemake.input.zarr)
write_store_to_zip(source_store, snakemake.output.ozx, version="0.5")
