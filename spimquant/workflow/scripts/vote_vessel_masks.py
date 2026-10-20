"""Majority vote of binary vessel masks (0/100) into one 0/100 consensus mask."""

import dask.array as da
import numpy as np

MASK_TRUE_VALUE = np.uint8(100)


def majority_vote(arrays):
    """Return MASK_TRUE_VALUE where more than half of `arrays` are > 0, else 0.

    With two arrays this is their intersection; with three, any two must agree.
    """
    votes = sum(da.asarray(a > 0).astype(np.uint8) for a in arrays)
    return da.where(votes * 2 > len(arrays), MASK_TRUE_VALUE, np.uint8(0)).astype(
        np.uint8
    )


if __name__ == "__main__":
    from dask.diagnostics import ProgressBar
    from dask_setup import get_dask_client
    from zarrnii import ZarrNii

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        masks = [ZarrNii.from_file(p) for p in snakemake.input.masks]
        shapes = {m.darr.shape for m in masks}
        if len(shapes) != 1:
            raise ValueError(f"Vessel masks to vote on differ in shape: {shapes}")

        out = masks[0]
        out.darr = majority_vote([m.darr for m in masks])

        with ProgressBar():
            # Match pyramid generation used by other mask-producing scripts.
            out.to_ome_zarr(
                snakemake.output.mask,
                match_scale_factors_from=snakemake.input.masks[0],
                **snakemake.config["zarrnii_out_kwargs"],
            )
