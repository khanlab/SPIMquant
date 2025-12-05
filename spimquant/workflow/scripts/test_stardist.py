if __name__ == "__main__":

    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=int(snakemake.threads / 2),  # or 32, depending on workload
        threads_per_worker=2,  # isolate GIL
        memory_limit="auto",  # or tune to your RAM
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)


    from zarrnii import ZarrNii

    import zarrnii
    znimg = ZarrNii.from_ome_zarr('/nfs/trident3/lightsheet/prado/mouse_app_lecanemab_ki3/bids/sub-AS134F1/micr/sub-AS134F1_sample-brain_acq-imaris4x_SPIM.ome.zarr')
    print(znimg)

    import nibabel as nib
    import numpy as np
    from stardist.models import StarDist3D
    from csbdeep.utils import normalize
    from stardist.data import test_image_nuclei_3d
    from pathlib import Path


    model = StarDist3D(None, name=snakemake.input.model_dir)

    img = test_image_nuclei_3d()



    labels,_ = model.predict_instances(normalize(img))


    out_dir_path = Path(snakemake.output.out_dir)
    out_dir_path.mkdir()

    nib.Nifti1Image(img,affine=np.eye(4,4)).to_filename(out_dir_path / 'img.nii')
    nib.Nifti1Image(labels,affine=np.eye(4,4)).to_filename(out_dir_path / 'labels.nii')

