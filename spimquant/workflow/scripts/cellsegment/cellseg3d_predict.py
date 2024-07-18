import os.path

import numpy as np
import cvpl_tools.persistence as persistence
import cvpl_tools.np_algs as np_algs
import dask.array as da
from dask.distributed import Client, progress
import time

import zarr
from ome_zarr.writer import write_image, write_labels, write_label_metadata
from dask.diagnostics import ProgressBar
from ome_zarr.scale import Scaler
from ome_zarr.io import parse_url
import argparse
import shutil


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        prog='cellseg3d_predict_ome_zarr',
        description='A program predicting on OME_ZARR lightsheet data yielding binary cell segmentation',
        epilog='This program will output a single label file of binary type with the same shape as the input, in OME_ZARR '
               'label-image format, but by default without copying over the accompanying input image.'
    )
    argparser.add_argument('--GPU', help='Use GPUs for computing',
                           action='store_true')
    argparser.add_argument('--USE_SYNTHETIC_DATASET', help='Generate synthetic dataset in the input path to '
                                                           'be used for testing, instead of predicting on real data',
                           action='store_true')
    argparser.add_argument('--COMPUTE_WIDTH', help='Width of input cube for CellSeg3D prediction',
                           type=int, default=64)
    argparser.add_argument('--NTHREAD', help='Number of threads per worker',
                           type=int)
    argparser.add_argument('--NWORKER', help='Number of workers',
                           type=int)
    argparser.add_argument('--COPY_INPUT', help='Copying input image over to the same directory',
                           action='store_true')
    argparser.add_argument('--CELLSEG3D_REPO_PATH', help='Repository path to CellSeg3D repo',
                           type=str)
    argparser.add_argument('--CELLSEG3D_CONFIG_PATH', help='Repository path to model config',
                           type=str)
    argparser.add_argument('--ZARR_PATH', help='Path to input OME_ZARR file',
                           type=str)
    argparser.add_argument('--OUT_ZARR_PATH', help='Path to output OME_ZARR file (to be created)',
                           type=str)
    argparser.add_argument('--CLASSIFY_CHANNEL', help='The channel to predict on',
                           type=int)
    argparser.add_argument('--TMP_PATH', help='where temporary files can be stored for program',
                           type=str)
    default_args = dict(snakemake.params.items())
    argparser.set_defaults(**default_args)
    args = argparser.parse_args()

    if args.COPY_INPUT and args.ZARR_PATH == args.OUT_ZARR_PATH:
        raise ValueError('ERROR: When args.COPY_INPUT=True, the output directory can not '
                         'be the same OME_ZARR file as input!')
    print('###argparse arguments received:\n', args.__dict__)  # print all arguments received
    if args.ZARR_PATH == '/path/to/input/ome/zarr':
        raise ValueError('Something went wrong in argparse. If your arguments are correctly received by the program')

    args.DEVICE = 'cuda' if args.GPU else 'cpu'

    args.CLASSIFY_CHANNEL = 0  # only predict on the first channel
    args.CHUNK_SIZE = (1, args.COMPUTE_WIDTH, 1024, 1024)


def write_da_as_ome_zarr_direct(zarr_group: zarr.Group, da_arr=None, lbl_arr=None, MAX_LAYER=3):
    if da_arr is not None:
        # assert the group is empty, since we are writing a new group
        for mem in zarr_group:
            raise ValueError('ZARR group is Non-empty, please remove the original ZARR before running the program to '
                             f'create synthetic data. ZARR group: {zarr_group}')

    scaler = Scaler(max_layer=MAX_LAYER, method='nearest')
    coordinate_transformations = []
    for layer in range(MAX_LAYER + 1):
        coordinate_transformations.append(
            [{'scale': [1., 1., (2 ** layer) * 1., (2 ** layer) * 1.],  # image-pyramids in XY only
              'type': 'scale'}])
    with ProgressBar():
        if da_arr is not None:
            write_image(image=da_arr,
                        group=zarr_group,
                        scaler=scaler,
                        coordinate_transformations=coordinate_transformations,
                        storage_options={'dimension_separator': '/'},
                        axes=['c', 'z', 'y', 'x'])

        # we could just use ['c', 'z', 'y', 'x'], however, napari ome zarr can't handle channel types but only space
        # type axes. So we need to fall back to manual definition, avoid 'c' which defaults to a channel type
        lbl_axes = [{'name': ch, 'type': 'space'} for ch in ['c', 'z', 'y', 'x']]
        if lbl_arr is not None:
            lbl_name = 'arrlbl'
            import numcodecs
            compressor = numcodecs.Blosc(cname='lz4', clevel=9, shuffle=numcodecs.Blosc.BITSHUFFLE)
            write_labels(labels=lbl_arr,
                         group=zarr_group,
                         scaler=scaler,
                         name=lbl_name,
                         coordinate_transformations=coordinate_transformations,
                         storage_options=dict(compressor=compressor),
                         axes=lbl_axes)
            # write_label_metadata(group=g,
            #                      name=f'/labels/{lbl_name}',
            #                      properties=properties)


def cache_image(arr: da.Array, location: str):
    store = parse_url(location, mode='w').store
    g = zarr.group(store)
    write_da_as_ome_zarr_direct(g, arr, None, 0)
    zarr_group = zarr.open(location, mode='r')
    zarr_subgroup = zarr_group['0']
    arr_read = da.from_zarr(zarr_subgroup)
    return arr_read


def write_da_as_ome_zarr(ome_zarr_path, da_arr=None, lbl_arr=None):
    store = parse_url(ome_zarr_path, mode='w').store
    g = zarr.group(store)
    if da_arr is not None:
        path1 = f'{args.TMP_PATH}/im1'
        if os.path.exists(path1):
            shutil.rmtree(path1)
        da_arr = cache_image(da_arr, path1)
    if lbl_arr is not None:
        path2 = f'{args.TMP_PATH}/im2'
        if os.path.exists(path2):
            shutil.rmtree(path2)
        lbl_arr = cache_image(lbl_arr, path2)
    write_da_as_ome_zarr_direct(g, da_arr, lbl_arr, 3)


def get_ome_zarr() -> da.Array:
    if args.USE_SYNTHETIC_DATASET:
        arr: da.Array = da.zeros((2, 255, 600, 600), dtype=np.uint16,
                                 chunks=(1, 1, 256, 256))

        def process_block(block, block_info=None):
            if block_info is not None:
                # calculate (global) indices array for each pixel
                block_slice = block_info[0]['array-location']
                indices = np.indices(block.shape)
                for dim in range(indices.shape[0]):
                    indices[dim] += block_slice[dim][0]
            else:
                print('block_info is None')
                return np.zeros(block.shape, dtype=np.uint16)
            # now, create balls in the block
            sq = np.zeros(block.shape)  # distance squared
            for dim in range(1, indices.shape[0]):  # every dim except channel dim which does not have distance
                sq += np.power(indices[dim], 2.) * .0002
            for dim in range(1, indices.shape[0]):  # every dim except channel dim which does not have distance
                indices[dim] %= 32
                sq += np.power(indices[dim] - 15.5, 2.)
            im = np.array(np.clip(1200. - sq * 15., 0., 1200.), dtype=np.uint16)
            return im

        arr = arr.map_blocks(process_block, dtype=np.uint16)
        write_da_as_ome_zarr(args.ZARR_PATH, da_arr=arr)

    zarr_group = zarr.open(args.ZARR_PATH, mode='r')
    zarr_subgroup = zarr_group['0']
    da_arr = da.from_zarr(zarr_subgroup)
    return da_arr


def classify_block(block, cmd_args, block_info=None):
    args = cmd_args

    import sys
    sys.path.append(args.CELLSEG3D_REPO_PATH)
    import napari_cellseg3d.create_model as create_model
    import napari_cellseg3d.predict as predict
    import torch
    import threading
    from dask.distributed import get_worker

    TID = threading.get_ident()
    if block_info is not None:
        # calculate (global) indices array for each pixel
        block_slice = block_info[0]['array-location']
        print(TID, 'start', time.time(), block_slice)
    else:
        return np.zeros(block.shape, dtype=np.bool_)
    worker = get_worker()
    if not hasattr(worker, 'worker_tup'):
        print(TID, 'creating tuple')
        config: dict = persistence.read_dict(args.CELLSEG3D_CONFIG_PATH)
        config['device'] = args.DEVICE
        cellseg3d_model = create_model.create_model(config)
        worker.worker_tup = (config, cellseg3d_model)
    else:
        print(TID, 'getting tuple')
    worker_tup = worker.worker_tup
    config, cellseg3d_model = worker_tup

    def get_mask(im):
        """
        Here im is a 3D (D * H * W) image;
        note on edges im may have a smaller width than args.COMPUTE_WIDTH in one of the dimensions, this can be
        accommodated by padding on that dimension
        """
        assert len(im.shape) == 3, f'Expected 3d image got shape {im.shape}'
        orig_shape = im.shape
        PADDING_FLAG = any([d < args.COMPUTE_WIDTH for d in orig_shape])
        if PADDING_FLAG:
            padded_dims = tuple(max(args.COMPUTE_WIDTH, d) for d in orig_shape)
            padded_im = np.zeros(padded_dims, dtype=im.dtype)
            padded_im[:im.shape[0], :im.shape[1], :im.shape[2]] = im
            im = padded_im

        im_slice_torch = torch.from_numpy(np.array(im * 1000., dtype=np.float32)[None, None])
        if args.GPU:
            CSIZE = (args.COMPUTE_WIDTH, im.shape[1], im.shape[2])
            seg = np.zeros((config['num_classes'], ) + im.shape, dtype=np.float32)
            for z in range(0, im_slice_torch.shape[2], CSIZE[0]):
                for y in range(0, im_slice_torch.shape[3], CSIZE[1]):
                    for x in range(0, im_slice_torch.shape[4], CSIZE[2]):
                        print(CSIZE, z, y, x, im_slice_torch.shape)
                        im_cc = im_slice_torch[:, :, z:z + CSIZE[0], y:y + CSIZE[1], x:x + CSIZE[2]].cuda()

                        seg_cc = predict.inference_on_torch_batch(config, im_cc,
                                                               [args.COMPUTE_WIDTH, args.COMPUTE_WIDTH, args.COMPUTE_WIDTH],
                                                               cellseg3d_model)
                        seg_cc = seg_cc.detach().cpu().numpy()[0]
                        seg[:, z:z + CSIZE[0], y:y + CSIZE[1], x:x + CSIZE[2]] = seg_cc
        else:
            seg = predict.inference_on_torch_batch(config, im_slice_torch,
                                                   [args.COMPUTE_WIDTH, args.COMPUTE_WIDTH, args.COMPUTE_WIDTH],
                                                   cellseg3d_model)
            seg = seg.numpy()[0]
        seg = seg[3]
        # seg = scipy.ndimage.gaussian_filter(seg, sigma=.62, mode='nearest')

        if PADDING_FLAG:
            seg_bin = seg[:orig_shape[0], :orig_shape[1], :orig_shape[2]] > .45
        else:
            seg_bin = seg > .45
        return seg_bin

    # counter = np_algs.get_counter(12, get_mask)
    #
    # def get_inst(block):
    #     inst = np.zeros(block.shape, dtype=np.uint16)
    #     for i in range(block.shape[0]):
    #         mask = counter.s1.mask(block[i])
    #         inst[i] = counter.s2.inst(mask)
    #     return inst

    mask = np.zeros(block.shape, dtype=np.bool_)
    for i in range(block.shape[0]):
        mask[i] = get_mask(block[i])
    print(TID, 'end', time.time())
    return mask


def main():
    import dask
    with dask.config.set({'temporary_directory': args.TMP_PATH}):
        if args.NWORKER > 1:
            import socket
            client = Client(f'{socket.gethostname()}:8786')
        else:
            client = Client(threads_per_worker=args.NTHREAD, n_workers=1)
        da_arr = get_ome_zarr()[args.CLASSIFY_CHANNEL:args.CLASSIFY_CHANNEL + 1]
        da_arr = da_arr.rechunk(chunks=args.CHUNK_SIZE)
        lbl_arr = da_arr.map_blocks(classify_block, dtype=np.bool_, cmd_args=args)

        copy_arr = da_arr if args.COPY_INPUT else None
        write_da_as_ome_zarr(args.OUT_ZARR_PATH, da_arr=copy_arr, lbl_arr=lbl_arr)  # if OUT_ZARR_PATH is different from input path, then only write a label


if __name__ == '__main__':
    main()
