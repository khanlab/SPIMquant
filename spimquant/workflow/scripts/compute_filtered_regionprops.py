import os
import shutil
import tempfile
import zipfile
from dask_setup import get_dask_client
from zarrnii import ZarrNii

if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        # 1. Create a secure, temporary directory in the system default tempdir
        with tempfile.TemporaryDirectory(suffix=".ome.zarr") as temp_dir:
            print(f"Extracting zip archive to temporary directory: {temp_dir}")

            # 2. Open and extract the entire input zip file safely
            with zipfile.ZipFile(snakemake.input.mask, "r") as zip_ref:
                zip_ref.extractall(temp_dir)

            # 3. Locate the extracted directory/file path inside the temp folder
            # OME-Zarr is usually a single top-level directory inside the zip.
            print("looking at extracted folder")
            extracted_contents = os.listdir(temp_dir)
            if not extracted_contents:
                raise ValueError("The input zip file is empty.")

            print("ready to read it")

            # Direct path to the extracted .zarr directory structure
            zarr_temp_path = os.path.join(temp_dir, extracted_contents[0])

            # 4. Point ZarrNii to the unzipped DirectoryStore instead of the ZipStore
            znimg = ZarrNii.from_file(
                temp_dir,
                level=0,  # input image is already downsampled to the wildcard level
            )

            znimg.compute_region_properties(
                output_path=snakemake.output.regionprops_parquet,
                region_filters=snakemake.params.region_filters,
                output_properties=snakemake.params.output_properties,
            )
