import os
import shutil
import glob
import numpy as np
import nibabel as nib
from PIL import Image
import scipy
from dataclasses import dataclass, field
from .array_key_dict import ArrayKeyDict

"""
fs.py
This is a file that provides methods to query and modify files in the file system
"""


def ensure_dir_exists(dir_path, remove_if_already_exists):
    if os.path.exists(dir_path):
        if remove_if_already_exists:
            shutil.rmtree(dir_path)
            os.mkdir(dir_path)
    else:
        os.mkdir(dir_path)


class ImFileType:
    FTYPE_UNKNOWN = -1
    FTYPE_JPX = 0
    FTYPE_JPG = 1
    FTYPE_PNG = 2
    FTYPE_NIB = 3
    FTYPE_TIF = 4
    FTYPE_NPY = 5
    FTYPE_NPZ = 6  # probabilities output by nnunet
    ext_to_ftype = {
        '.jpx': FTYPE_JPX,
        '.jpg': FTYPE_JPG,
        '.png': FTYPE_PNG,
        '.nii.gz': FTYPE_NIB,
        '.tif': FTYPE_TIF,
        '.npy': FTYPE_NPY,
        '.npz': FTYPE_NPZ,
    }
    ftype_to_ext = list(ext_to_ftype.keys())
    all_ftypes = list(ext_to_ftype.values())  # except unknown

    FORMAT_UINT8 = 0
    FORMAT_UINT16 = 1
    FORMAT_INT32 = 2
    FORMAT_FLOAT32 = 3
    format_to_np_dtype = [np.uint8, np.uint16, np.int32, np.float32]

    @classmethod
    def ftype_from_im_path(cls, im_path: str) -> int:
        """
        Infer the ftype from image path;
        The path can contain wildcards like "*"
        When the ftype can't be inferred, FTYPE_UNKNOWN will be returned
        """
        prefix, ext = os.path.splitext(im_path)
        ext = ext.lower()
        if ext == '.gz':
            prefix_prefix, prefix_ext = os.path.splitext(prefix)
            if prefix_ext.lower() == '.nii':
                return ImFileType.FTYPE_NIB
            else:
                return ImFileType.FTYPE_UNKNOWN
        elif ext in ImFileType.ext_to_ftype:
            return ImFileType.ext_to_ftype[ext]
        else:
            return ImFileType.FTYPE_UNKNOWN

    @classmethod
    def get_imid_and_zstack(cls, im_path: str, is_stack=False, keep_dirpath=False) -> tuple[str, str]:
        """
        Get the path stem as the unique image id from the path of the file;
        Second return value is the zstack number as a string
        file extensions will be stripped;
        If the file is part of a stacked image, then the z-stack numbers will be stripped from stem as will
        """
        ftype = ImFileType.ftype_from_im_path(im_path)
        if ftype == ImFileType.FTYPE_UNKNOWN:
            raise ValueError(f'Attempt to get the imid of an image file from path "{im_path}" with unknown file type!')
        elif ftype == ImFileType.FTYPE_NIB:
            imid = im_path[:-len('.nii.gz')]
        else:
            imid = os.path.splitext(im_path)[0]

        if not keep_dirpath:
            imid = os.path.split(imid)[1]

        if is_stack:
            i = -1
            while len(imid) >= -i and imid[i].isdigit():
                i -= 1
            return imid[:i + 1], imid[i + 1:]
        else:
            return imid, ''

    @classmethod
    def get_imid(cls, im_path: str, is_stack=False, keep_dirpath=False) -> str:
        return ImFileType.get_imid_and_zstack(im_path, is_stack, keep_dirpath)[0]

    @classmethod
    def load_im_as_np(cls, im_format, path, true_im_ndim: int, stack_axis=None,
                      concat_instead=False, allow_pickle=False) -> tuple:
        """
        Load individual file as a np array, or a stack of images as a np array if a list of path is given
        Returns a tuple (im, im_ftype, im_meta) where im_meta may be needed in order to save the im back to a file
        for certain file types
        """
        if isinstance(path, str):
            ftype = ImFileType.ftype_from_im_path(path)
            assert stack_axis is None, 'ERROR: Passing stack_axis in when loading single file'
            assert ftype != ImFileType.FTYPE_UNKNOWN, 'ERROR: Attempt to load an unknown image type!'
            if ftype == ImFileType.FTYPE_NPZ:
                assert im_format == ImFileType.FORMAT_FLOAT32
                im = np.load(path, allow_pickle=allow_pickle)['probabilities']
                im = np.transpose(im, (0, 3, 2, 1))
                im_meta = [ftype, im.shape, None, path, None]
            elif ftype == ImFileType.FTYPE_NPY:
                im = np.load(path, allow_pickle=allow_pickle)
                if allow_pickle:
                    im = im.item()['masks']
                im_meta = [ftype, im.shape, None, path, None]
            elif ftype == ImFileType.FTYPE_NIB:
                nii_im = nib.load(path)
                im = np.array(nii_im.get_fdata(), dtype=ImFileType.format_to_np_dtype[im_format])
                im_meta = [ftype, im.shape, None, path, nii_im.affine]
            else:
                im = np.array(Image.open(path), dtype=ImFileType.format_to_np_dtype[im_format])
                im_meta = [ftype, im.shape, None, path, None]

            assert im.dtype == ImFileType.format_to_np_dtype[im_format], \
                f'ERROR: Expect file type {ImFileType.format_to_np_dtype[im_format]} but got an image of dtype {im.dtype} at {path}'
            if len(im.shape) > true_im_ndim:
                axes = range(true_im_ndim, len(im.shape))
                for i in axes:
                    assert im.shape[i] == 1, f'Extra axes to true_im_ndim must be of width 1. true_im_ndim={true_im_ndim}, im.shape={im.shape}'
                im = im.squeeze(axis=tuple(axes))
            assert len(im.shape) == true_im_ndim, (f'ERROR: Expected true dimension {true_im_ndim}, but found '
                                                  f'image of shape {im.shape}!')
            assert type(im_meta) is list, 'ERROR: im_meta is not properly initialized!'
            return im, im_meta
        else:
            assert isinstance(path, list)
            im_arr = []
            im_meta = None
            for filepath in path:
                if concat_instead:
                    true_slice_dim = true_im_ndim
                else:
                    true_slice_dim = true_im_ndim - 1
                tup = ImFileType.load_im_as_np(im_format, filepath, true_slice_dim)
                im_arr.append(tup[0])
                cur_im_meta = tup[1]
                if im_meta is None:
                    im_meta = list(cur_im_meta)  # shallow copy
                    ftype_meta = im_meta[-1]
                    im_meta[-1] = [ftype_meta]
                else:
                    assert im_meta[0] == cur_im_meta[0]  # same file type
                    assert im_meta[1] == cur_im_meta[1]  # each image has the same shape
                    im_meta[-1].append(cur_im_meta[-1])
            if stack_axis is None:
                raise ValueError(f"Please specify stack_axis!")
                # stack_axis = true_im_ndim - 1  # default stack on the last axis
            if concat_instead:  # concat on existing axis
                assert 0 <= stack_axis < true_im_ndim - 1, f'Attempt to concat on axis {stack_axis} with true dim={true_im_ndim}'
                im = np.concatenate(im_arr, axis=stack_axis)
            else:  # stack on new axis
                assert 0 <= stack_axis < true_im_ndim, f'Attempt to concat on axis {stack_axis} with true dim={true_im_ndim}'
                im = np.stack(im_arr, axis=stack_axis)

            im_meta[2] = stack_axis
            im_meta[3] = path

            return im, im_meta

    @classmethod
    def save_np_as_im(cls, target_path: str, im: np.ndarray, ftype=None, write_axis_order=None, stack_axis=None,
                      concat_instead=False, ftype_meta=None, allow_pickle=False, im_format=None):
        """
        If ftype is None (default), it will be inferred from path
        """
        true_im_ndim = len(im.shape)
        if type(target_path) is str:
            assert not os.path.exists(target_path), f'{target_path} exists!'
            if ftype is None:
                ftype = ImFileType.ftype_from_im_path(target_path)
            if im_format is not None:
                assert im.dtype == ImFileType.format_to_np_dtype[im_format]
            assert ftype != ImFileType.FTYPE_UNKNOWN, 'ERROR: Attempt to save an unknown image type!'
            if ftype == ImFileType.FTYPE_NPZ:
                assert im.dtype == np.float32
                im = np.transpose(im, (0, 3, 2, 1))
                np.savez(target_path, probabilities=im, allow_pickle=allow_pickle)
            elif ftype == ImFileType.FTYPE_NPY:
                if allow_pickle:
                    im = {'masks': im}
                np.save(target_path, im, allow_pickle=allow_pickle)
            elif ftype == ImFileType.FTYPE_NIB:
                assert ftype_meta.shape == (4, 4), f'Affine should have shape 4, 4, but get {ftype_meta[-1].shape}'
                nii_im = nib.Nifti1Image(im, ftype_meta)
                nib.save(nii_im, target_path)
            else:
                Image.fromarray(im).save(target_path)
        else:
            if stack_axis is None:
                stack_axis = true_im_ndim - 1
            for i in range(im.shape[stack_axis]):
                if concat_instead:
                    im_slice = im.take(indices=range(i, i + 1), axis=stack_axis)
                else:
                    im_slice = im.take(indices=i, axis=stack_axis)
                target = target_path[i]
                if ftype_meta is None:
                    cur_ftype_meta = None
                else:
                    cur_ftype_meta = ftype_meta[i]
                # recursive call!
                ImFileType.save_np_as_im(
                    target_path=target,
                    im=im_slice,
                    ftype=ftype,
                    write_axis_order=write_axis_order,
                    stack_axis=stack_axis,
                    concat_instead=False,
                    ftype_meta=cur_ftype_meta,
                    allow_pickle=allow_pickle)


@dataclass
class ImReadSetting:
    """
    im_format (int) - specifies the data storage layout e.g. ImFileType.FORMAT_UINT8
    stack_axis (int) - if not None, combine several files with different zstacks to the same imid; specify the axis to stack on
    concat_instead (bool) - if True, the stacking are done by np.concatenate() which preserves im dimensions; over the same axis
    true_im_ndim (int) - # of dimensions of a single image; includes number of channels (for 2d BGR image this would be 3)
    allow_pickle (bool) - set allow_pickle for np.load and np.save
    imid_keep_dirpath (bool) - keeping the directory path when extracting imid from paths; default to False
    """
    im_format: int = ImFileType.FORMAT_UINT8
    stack_axis: int = None
    concat_instead: bool = False
    true_im_ndim: int = 3
    allow_pickle: bool = False
    imid_keep_dirpath: bool = False


@dataclass
class ImWriteSetting:
    """
    im_format (int) - specifies the data storage layout e.g. ImFileType.FORMAT_UINT8
    stack_axis (int) - if not None, combine several files with different zstacks to the same imid; specify the axis to stack on
    concat_instead (bool) - if True, the stacking are done by np.concatenate() which preserves im dimensions; over the same axis
    allow_pickle (bool) - set allow_pickle for np.load and np.save
    ftype (int) - file type to write e.g. ImFileType.FTYPE_NPZ
    """
    im_format: int = ImFileType.FORMAT_UINT8
    stack_axis: int = None
    concat_instead: bool = False
    allow_pickle: bool = False
    ftype: int = None


class ImIO:
    @classmethod
    def read_single_image(cls, read_setting: ImReadSetting, path):
        im, im_meta = ImFileType.load_im_as_np(read_setting.im_format, path,
                                               true_im_ndim=read_setting.true_im_ndim,
                                               stack_axis=read_setting.stack_axis,
                                               concat_instead=read_setting.concat_instead,
                                               allow_pickle=read_setting.allow_pickle)
        return im, im_meta

    @classmethod
    def write_single_image(cls, write_setting: ImWriteSetting, path, im, im_meta):
        ImFileType.save_np_as_im(path, im=im, ftype=write_setting.ftype,
                                 write_axis_order=write_setting.write_axis_order,
                                 stack_axis=write_setting.stack_axis,
                                 concat_instead=write_setting.concat_instead,
                                 ftype_meta=im_meta[-1],
                                 allow_pickle=write_setting.allow_pickle)

    @classmethod
    def read_filenames(cls, read_setting: ImReadSetting, pattern: str = None, path=None):
        """
        The keys (imids) of the ArrayKeyDict returned by this function will be sorted alphabetically
        Parameters
            path (list) - optional and passed in as a list of all paths if pattern is None
        Return
            imid_to_path (Dict) - maps from imids to either a filename or a list of filenames (if image is stacked);
                The keys of this dictionary is a sorted list of image imids, act as unique ids of images, sorted alphabetically
        """
        imid_to_path = ImIO._get_imid_to_path_mapping(read_setting, path, pattern)
        imid_to_path.reorder_keys(sorted(imid_to_path.keys()))
        return imid_to_path

    @classmethod
    def read_image(cls, read_setting: ImReadSetting, paths: list):
        """
        Parameter
            paths (list) - can be nested in case of stacked images
        Return
            an iterator whose elements are (im, im_meta), im is the np array of the image and im_meta contains much
            additional metadata returned by ImFileType.load_im_as_np()
        """
        assert type(paths) is not ArrayKeyDict, f'Attempt to pass an ArrayKeyDict to read_image()!'
        for path in paths:
            im, im_meta = ImIO.read_single_image(read_setting, path)
            yield im, im_meta

    @classmethod
    def write_image(cls, write_setting: ImWriteSetting, paths: list, im_iter: iter):
        """
        Parameter
            paths (list) - can be nested in case of stacked images
            im_iter (iter) - An iterator returning tuples (im, im_meta), the order of these images should be the same as
                the order of their paths in the paths variable
        """
        for path in paths:
            im, im_meta = next(im_iter)
            ImIO.write_single_image(write_setting, path, im, im_meta)

    @classmethod
    def _get_imid_to_path_mapping(cls, read_setting: ImReadSetting,
                                  paths: list[str] = None, pattern: str = None):
        if paths is None:
            paths = glob.glob(pattern)
        elif pattern is not None:
            print('WARNING: both paths and pattern are not None! Is this intended?')

        imid_to_path = ArrayKeyDict()
        for path in paths:
            ftype = ImFileType.ftype_from_im_path(path)
            if ftype == ImFileType.FTYPE_UNKNOWN:
                continue
            imid, zstack_str = ImFileType.get_imid_and_zstack(path, is_stack=read_setting.stack_axis is not None,
                                                              keep_dirpath=read_setting.imid_keep_dirpath)
            if zstack_str == '':
                imid_to_path[imid] = path
            else:
                tup = (zstack_str, path)
                if imid not in imid_to_path:
                    imid_to_path[imid] = [tup]
                else:
                    imid_to_path[imid].append(tup)
        for k in list(imid_to_path.keys()):
            v = imid_to_path[k]
            if type(v) is list:
                imid_to_path[k] = [tup[1] for tup in sorted(v)]
        return imid_to_path


def im_resize(target_size: tuple, im: np.ndarray, order: int = 0) -> np.ndarray:
    assert len(im.shape) == len(target_size), f'Image has shape {im.shape} but intend to resize ' \
                                              f'to incompatible shape {target_size}'
    zoom_factors = (target_size[i] / im.shape[i] for i in range(len(target_size)))
    resized_im = scipy.ndimage.zoom(im, zoom_factors, order=order)  # order=0 -> use nearest neighbor
    return resized_im


def test():
    assert ImFileType.ftype_from_im_path('./some_root/hello.nii.gz') == ImFileType.FTYPE_NIB
    assert ImFileType.ftype_from_im_path('./some_root/hello_0001.JPG') == ImFileType.FTYPE_JPG
    assert ImFileType.ftype_from_im_path('./a.py') == ImFileType.FTYPE_UNKNOWN
    assert ImFileType.get_imid('./some_root/hello.png') == 'hello'
    assert (ImFileType.get_imid('./some_root/nibfile0001.nii.gz', is_stack=True, keep_dirpath=True) ==
            './some_root/nibfile')

    print('All tests successfully completed')
