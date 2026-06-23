"""Rechunk SPIM OME-Zarr to work directory for optimized processing.

Copies all channels and scale levels from the input OME-Zarr to the work
directory with configurable chunking.  The output is always an ome.zarr
directory on the local work disk.

Zarr groups are traversed recursively: data arrays are rechunked using dask,
and all metadata attributes are preserved verbatim.
"""

import zarr
import dask.array as da
from zarr.storage import LocalStore

from dask_setup import get_dask_client


def _open_src_store(path):
    """Open a zarr-compatible store for various input formats."""
    if path.endswith(".ozx") or path.endswith(".ome.zarr.zip"):
        return zarr.storage.ZipStore(path)
    return LocalStore(path)


def _copy_rechunked(src_grp, dst_grp, chunks):
    """Recursively copy a zarr group, rechunking data arrays.

    Attributes are preserved on all groups and arrays.  Arrays are rechunked to
    the requested size (clamped to the array shape so no empty chunks exist at
    boundaries).  If the chunks tuple has fewer dimensions than the array, the
    array's existing chunk sizes are used for the remaining dimensions.
    """
    dst_grp.attrs.update(dict(src_grp.attrs))
    for key in src_grp:
        item = src_grp[key]
        if isinstance(item, zarr.Array):
            n_dims = len(item.shape)
            # Extend chunks to match array dimensionality using existing chunk
            # sizes for any extra dimensions not covered by the config value.
            padded_chunks = list(chunks) + list(item.chunks[len(chunks) :])
            actual_chunks = tuple(
                min(c, s) for c, s in zip(padded_chunks[:n_dims], item.shape)
            )
            da_arr = da.from_zarr(item)
            dst_arr = dst_grp.create_array(
                key,
                shape=item.shape,
                chunks=actual_chunks,
                dtype=item.dtype,
                overwrite=True,
            )
            da.store(da_arr.rechunk(actual_chunks), dst_arr)
            dst_arr.attrs.update(dict(item.attrs))
        elif isinstance(item, zarr.Group):
            _copy_rechunked(item, dst_grp.require_group(key), chunks)


if __name__ == "__main__":
    input_path = str(snakemake.input.spim)
    output_path = str(snakemake.output.zarr)
    chunks = tuple(snakemake.params.chunks)

    if input_path.endswith(".ims"):
        raise ValueError(
            "--use-work-zarr does not support Imaris (.ims) format inputs. "
            "Convert to OME-Zarr format before using this option."
        )

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        src_store = _open_src_store(input_path)
        src_group = zarr.open_group(src_store, mode="r")

        dst_store = LocalStore(output_path)
        # Use zarr_format=2 to match the OME-NGFF v0.4/v0.5 format expected by ZarrNii.
        dst_group = zarr.open_group(dst_store, mode="w", zarr_format=2)

        try:
            _copy_rechunked(src_group, dst_group, chunks)
            zarr.consolidate_metadata(dst_store)
        finally:
            # Close zip stores explicitly to flush file handles; directory stores
            # do not require explicit closing but we call close() defensively.
            for store in (src_store, dst_store):
                if hasattr(store, "close"):
                    store.close()
