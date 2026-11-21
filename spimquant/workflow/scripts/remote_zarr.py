"""OME-Zarr outputs kept in object storage (MIND), with a local pointer file for Snakemake.

With --zarr_store_root, every OME-Zarr output of the workflow (``bids_oz_out``) is a small local
file ending in ``.zarr.url`` whose content is the URL of the image, written under that root, e.g.
``mind://scratch/study/derivatives/spimquant/<run>/sub-01/.../sub-01_..._SPIM.<id>.ome.zarr``.
Snakemake tracks the pointer files; the images themselves are read and written by URL.

MIND refuses to overwrite or delete objects, so each image written gets a key of its own (``<id>``
is new each time): a rerun of a job, or a retry of a failed one, never writes over an earlier
attempt. Earlier attempts and ``temp()`` outputs stay in the store, under the run's prefix.

Importing this module makes ``ZarrNii.from_file``, ``ZarrNii.from_ome_zarr`` and
``ZarrNii.to_ome_zarr`` understand pointer files, so scripts keep passing ``snakemake.input`` and
``snakemake.output`` paths. Other paths are handled as before.
"""

import os
import posixpath
import re
import uuid
from contextlib import contextmanager

from zarr.storage import WrapperStore

POINTER_SUFFIX = ".zarr.url"
# set by the workflow (Snakefile) when --zarr_store_root is given
ROOT_ENV = "SPIMQUANT_ZARR_STORE_ROOT"


def is_pointer(path) -> bool:
    return isinstance(path, (str, os.PathLike)) and os.fspath(path).endswith(
        POINTER_SUFFIX
    )


def resolve(path):
    """The URL a pointer file holds, or path itself if it is not a pointer file."""
    if not is_pointer(path):
        return path
    with open(path) as f:
        return f.read().strip()


def new_url(pointer_path) -> str:
    """A new URL for the image a pointer file will point to, under the run's root."""
    root = os.environ.get(ROOT_ENV)
    if not root:
        raise RuntimeError(
            f"{pointer_path} is a pointer file, but {ROOT_ENV} is not set"
        )
    rel = os.path.relpath(os.path.abspath(pointer_path)).replace(os.sep, "/")
    if rel.startswith("../"):
        raise ValueError(
            f"pointer file {pointer_path} is outside the working directory"
        )
    # MIND's proxy rejects keys with "+" or spaces (e.g. desc-gmm+n3k1); the pointer
    # file records the URL, so the key need not match the local name exactly
    stem = re.sub(r"[^A-Za-z0-9._/-]", "-", rel[: -len(POINTER_SUFFIX)])
    return posixpath.join(root.rstrip("/"), f"{stem}.{uuid.uuid4().hex[:8]}.ome.zarr")


@contextmanager
def open_for_writing(url):
    """A zarr group to write an OME-Zarr image into at url, completed on leaving the block."""
    import zarr

    if url.startswith("mind://"):
        from mind_client.base import CachedToken
        from mind_client.fs import DEFAULT_URL, key_of
        from mind_client.zarr_store import zarr_store

        server = os.environ.get("MIND_URL") or DEFAULT_URL
        # metadata is held back and uploaded once on leaving the block, since MIND refuses
        # overwrites; the token is reread from MIND_TOKEN_FILE, so dask workers can write too
        with zarr_store(
            server, key_of(url), CachedToken(os.environ.get("MIND_TOKEN_FILE"))
        ) as store:
            yield zarr.open_group(NoDeleteStore(store), mode="a")
    else:
        yield zarr.open_group(url, mode="w")


class NoDeleteStore(WrapperStore):
    """Wraps a store that refuses deletes (MIND): deleting an object that is not there is a no-op.

    zarr deletes a chunk's object when the chunk holds only the fill value, in case an earlier
    write left one; each image here is written under a new key, so there is never one.
    """

    def _with_store(self, store):
        return type(self)(store)

    async def delete(self, key):
        if await self._store.exists(key):
            raise PermissionError(f"cannot delete {key}: the store refuses deletes")


def _install():
    from zarrnii import ZarrNii

    if getattr(ZarrNii, "_spimquant_remote_zarr", False):
        return
    from_file = ZarrNii.from_file.__func__
    from_ome_zarr = ZarrNii.from_ome_zarr.__func__
    to_ome_zarr = ZarrNii.to_ome_zarr

    def from_file_resolved(cls, path, **kwargs):
        if is_pointer(path):
            return from_ome_zarr(cls, resolve(path), **kwargs)
        return from_file(cls, path, **kwargs)

    def from_ome_zarr_resolved(cls, store_or_path, *args, **kwargs):
        return from_ome_zarr(cls, resolve(store_or_path), *args, **kwargs)

    def to_ome_zarr_pointer(self, store_or_path, *args, **kwargs):
        if "match_scale_factors_from" in kwargs:
            kwargs["match_scale_factors_from"] = resolve(
                kwargs["match_scale_factors_from"]
            )
        if not is_pointer(store_or_path):
            return to_ome_zarr(self, store_or_path, *args, **kwargs)
        url = new_url(store_or_path)
        with open_for_writing(url) as group:
            result = to_ome_zarr(self, group, *args, **kwargs)
        # the pointer is written last, once the image is complete
        tmp = f"{store_or_path}.tmp"
        with open(tmp, "w") as f:
            f.write(url + "\n")
        os.replace(tmp, store_or_path)
        print(f"wrote {url}", flush=True)
        return result

    ZarrNii.from_file = classmethod(from_file_resolved)
    ZarrNii.from_ome_zarr = classmethod(from_ome_zarr_resolved)
    ZarrNii.to_ome_zarr = to_ome_zarr_pointer
    ZarrNii._spimquant_remote_zarr = True


_install()
