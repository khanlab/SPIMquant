"""
This file finds, globally in a 3d brightness image, all the contours' size and location.

Inputs to this file:
1. Path to an ome-zarr 3d single channel brightness image
2. Path to a tiffile 3d binary mask that mask off regions of the image in 1) with false positives
3. A scaling factor tuple, when multiplied from left, transforms pixel location of 2) to 1)
4. An integer (starting from 0) specifies the index of input channel to use

Output of this file is a NDBlock list of centroids saved in local file, which can be loaded
and read as a global list of all contours' size and location.
"""


if __name__ == '__main__':
    import cvpl_tools.im.process.qsetup as qsetup
    # IMPORT YOUR LIBRARIES HERE
    import cvpl_tools.im.seg_process as seg_process
    import cvpl_tools.im.process.bs_to_os as sp_bs_to_os
    import cvpl_tools.im.process.os_to_cc as sp_os_to_cc
    import cvpl_tools.im.process.any_to_any as sp_any_to_any
    import cvpl_tools.im.algs.dask_resize as im_resize
    import cvpl_tools.ome_zarr.io as ome_io
    import cvpl_tools.im.ndblock as ndblock
    import dask.array as da
    import numcodecs
    import tifffile
    import shutil

    class Pipeline(seg_process.SegProcess):
        def __init__(self):
            super().__init__()
            self.in_to_bs = seg_process.SimpleThreshold(.45)
            self.bs_to_os = sp_bs_to_os.DirectBSToOS(is_global=True)
            self.os_to_cc = sp_os_to_cc.CountOSBySize(
                size_threshold=200.,
                volume_weight=5.15e-3,
                border_params=(3., -.5, 2.3),
                min_size=8,
                reduce=False,
                is_global=True
            )

        def forward(self, im, cptr, viewer_args: dict = None):
            cdir = cptr.subdir()
            bs = self.in_to_bs.forward(im, cptr=cdir.cache(cid='in_to_bs'), viewer_args=viewer_args)
            os = self.bs_to_os.forward(bs, cptr=cdir.cache(cid='bs_to_os'), viewer_args=viewer_args)
            cc = self.os_to_cc.forward(os, cptr=cdir.cache(cid='os_to_cc'), viewer_args=viewer_args)
            return cc

    TMP_PATH = snakemake.params.tmp_path
    with qsetup.PLComponents(TMP_PATH, 'CacheDirectoryNegMaskedCounting',
                             client_args=dict(threads_per_worker=12, n_workers=1),
                             viewer_args=dict(use_viewer=False)) as plc:
        # DO DASK COMPUTATION, AND SHOW RESULTS IN plc.viewer
        src_im = ome_io.load_dask_array_from_path(snakemake.params.zarr, mode='r', level=0)
        pipeline = Pipeline()
        src_im = da.clip(src_im / 1000, 0., 1.)
        assert src_im.ndim == 3
        print(f'Saving results in {plc.cache_root.abs_path}')
        print(f'Computing centroids size and location. Masking the image, imshape={src_im.shape}.')
        src_im = src_im.rechunk(chunks=(128, 256, 256))
        storage_options = dict(
            dimension_separator='/',
            preferred_chunksize=None,  # don't re-chunk when saving and loading
            multiscale=0,
            compressor=numcodecs.Blosc(cname='lz4', clevel=9, shuffle=numcodecs.Blosc.BITSHUFFLE)
        )
        viewer_args = dict(
            viewer=None,
            display_points=False,
            display_checkerboard=False,
            client=plc.dask_client,
            storage_options=storage_options
        )

        def compute_masking():
            neg_mask = da.from_array(tifffile.imread(snakemake.input.neg_mask), chunks=(64, 64, 64))
            neg_mask = im_resize.upsample_pad_crop_fit(
                src_arr=neg_mask,
                tgt_arr=src_im,
                cptr=plc.cache_root.cache('neg_mask_upsampling'),
                viewer_args=viewer_args | dict(is_label=True),
            )
            return src_im * (1 - neg_mask)
        plc.cache_root.cache_im(compute_masking, cid='masked_src_im', viewer_args=viewer_args)
        cc = pipeline.forward(src_im, plc.cache_root.cache(cid='global_label'), viewer_args=viewer_args)
        object_scores = cc.reduce(force_numpy=True)
        print(f'Total number of object score (not n_object nor volume) estimated: {object_scores.sum().item()}')

        tgt_folder = snakemake.output.found_lc
        shutil.move(f'{plc.cache_root.abs_path}/dir_cache_global_label/dir_cache_os_to_cc/dir_cache_os_to_lc/file_cache_lc_ndblock',
                    f'{tgt_folder}')
        lc_snapshot = ndblock.NDBlock.load(tgt_folder).reduce(force_numpy=True)[:20]
        print('Snapshot of list of centroids\n', lc_snapshot)
