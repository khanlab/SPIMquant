import skimage
import skimage.morphology as morph
import scipy.ndimage as ndimage
import numpy as np


def redistribute_weights(arr, ch):
    weighted = arr + .001
    arr[:] += arr[ch, None] * (weighted / (weighted.sum(axis=0) - weighted[ch])[None])
    arr[ch] = 0


def instance_to_binary(seg_inst):
    D, H, W = seg_inst.shape
    seg_bin = np.zeros((D - 1, H - 1, W - 1), dtype=np.uint8)
    chunk = seg_inst[:D - 1, :H - 1, :W - 1]
    for d in range(3):
        mins = [0, 0, 0]
        maxs = [D - 1, H - 1, W - 1]
        mins[d] += 1
        maxs[d] += 1
        chunk_other = seg_inst[mins[0]:maxs[0], mins[1]:maxs[1], mins[2]:maxs[2]]
        edge_detected = (chunk == chunk_other) + (chunk == 0) + (chunk_other == 0) == 0
        seg_bin = np.logical_or(seg_bin, edge_detected)
    seg_inst = np.copy(seg_inst)
    seg_inst[:D - 1, :H - 1, :W - 1] = np.logical_and(seg_inst[:D - 1, :H - 1, :W - 1], ~seg_bin)
    return seg_inst


def watershed(seg_bin, dist_thres=1., remove_smaller_than=None):
    """
    Run Watershed algorithm to perform instance segmentation. The result is a index labeled int64 mask
    :param dist_thres:
    :param seg_bin:
    :param remove_smaller_than:
    :return:
    """
    # reference: https://docs.opencv.org/4.x/d3/db4/tutorial_py_watershed.html
    fp_width = 2
    fp = [(np.ones((fp_width, 1, 1)), 1), (np.ones((1, fp_width, 1)), 1), (np.ones((1, 1, fp_width)), 1)]
    # sure_bg = morph.binary_dilation(seg_bin, fp)
    sure_bg = seg_bin
    # sure_fg = morph.binary_erosion(seg_bin, fp)
    dist_transform = ndimage.distance_transform_edt(seg_bin)
    sure_fg = dist_transform >= dist_thres  # doesn't work well: pixel around dist_thres=1 either all get cut or all retain
    unknown = sure_bg ^ sure_fg
    lbl_im = morph.label(sure_fg, connectivity=1)

    # so that sure_bg is 1 and unknown region is 0
    lbl_im += 1
    lbl_im[unknown == 1] = 0
    result = skimage.segmentation.watershed(-dist_transform, lbl_im, connectivity=1)
    result -= result > 0  # we don't need to mark sure_bg as 1, make all marker id smaller by 1

    if remove_smaller_than is not None:
        # watershed again over the small object regions
        small_mask, big_mask = split_labeled_objects(result, remove_smaller_than, connectivity=1)
        lbl_im[small_mask] = 0
        result = skimage.segmentation.watershed(-dist_transform, lbl_im, connectivity=1)

        # at this point, there may be unfilled space, which are rare cases where lots of small objects make up a
        # larger connected part; these space should be cells but the individual objects do not meet size requirement
        unfilled = morph.label(result == 0, connectivity=1)
        result += (unfilled != 0) * (unfilled + result.max())

        result -= result > 0  # we don't need to mark sure_bg as 1, make all marker id smaller by 1

    return result


def split_labeled_objects(lbl_im, size_thres, connectivity=1):
    component_sizes = np.bincount(lbl_im.ravel())
    small_inds = component_sizes < size_thres
    small_inds[0] = False
    small_mask = small_inds[lbl_im]
    big_inds = component_sizes >= size_thres
    big_inds[0] = False
    big_mask = big_inds[lbl_im]
    return small_mask, big_mask


def split_objects(seg, size_thres, connectivity=1):
    lbl_im = morph.label(seg, connectivity=connectivity)
    return split_labeled_objects(lbl_im, size_thres, connectivity)


def round_object_detection(seg, size_thres, dist_thres, rst, size_thres2, dist_thres2, rst2):
    """
    detect round objects in a 3d image
    :param seg: the binary mask where we want to detect on
    :param size_thres: the size threshold below which size, contours are all kept as is
    :param dist_thres: threshold for seeding in the watershed algorithm
    :return: same shape image labeled 0, 1, ..., n, where 0 is background and 1...n are detected objects
    """
    # objects too small cannot be connected. We use this property to first find small objects that must be single cells
    small_mask, big_mask = split_objects(seg, size_thres, connectivity=1)
    big_mask, big_mask2 = split_objects(big_mask, size_thres2, connectivity=1)

    # big objects may be a single large cell or overlapping small cells, run watershed to separate overlapping cells
    lbl_im = watershed(big_mask, dist_thres=dist_thres, remove_smaller_than=rst)
    lbl_im2 = watershed(big_mask2, dist_thres=dist_thres2, remove_smaller_than=rst2)

    # finally, combine the two set of cells
    small_labeled = morph.label(small_mask, connectivity=1)
    lbl_im += (small_labeled != 0) * (small_labeled + lbl_im.max())
    lbl_im += (lbl_im2 != 0) * (lbl_im2 + lbl_im.max())
    # lbl_im = (lbl_im2 > 0) * 1 + (lbl_im > 0) * 2 + small_mask * 3

    return lbl_im
