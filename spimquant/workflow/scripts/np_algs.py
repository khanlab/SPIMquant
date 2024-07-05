import enum
import skimage
import skimage.morphology as morph
import scipy.ndimage as ndimage
import numpy as np
from abc import ABC, abstractmethod
import skimage.draw as draw
import matplotlib.pyplot as plt


"""-------------------------Part 0: The interfaces----------------------------------------"""


class CountingMethod(enum.Enum):
    THRES_COUNT = 1
    THRES_BYSIZE = 2
    THRES_VOLUME = 3
    THRES_WATERSHED_COUNT = 4
    THRES_WATERSHED_BYSIZE = 5
    THRES_WATERSHED_VOLUME = 6
    CELLSEG3D_COUNT = 7
    CELLSEG3D_BYSIZE = 8
    CELLSEG3D_VOLUME = 9
    CELLSEG3D_WATERSHED_COUNT = 10
    CELLSEG3D_WATERSHED_BYSIZE = 11
    CELLSEG3D_WATERSHED_VOLUME = 12
    SUM_INTENSITY = 13
    BLOBDOG = 14


counting_method_inverse_dict = {item.value: item.name for item in CountingMethod}


def get_counter(ty: int, seg):
    if ty <= 12:
        if ty <= 6:
            s1 = Mask_Thres(.5)
        else:
            s1 = Mask_CellSeg3D(seg)
        VOXEL_PER_CELL = 60
        if ty % 3 == 0:
            ncell_from_inst = NCellFromInst_BySize(0., VOXEL_PER_CELL)
        elif ty % 3 == 1:
            ncell_from_inst = NCellFromInst_BySize(50., VOXEL_PER_CELL)
        else:
            ncell_from_inst = NCellFromInst_BySize(1e10, VOXEL_PER_CELL)

        use_watershed = ty % 6 >= 3
        if use_watershed:
            # s2 = CountFromMask_Watershed(ncell_from_inst,
            #                              size_thres,
            #                              dist_thres,
            #                              rst,
            #                              size_thres2,
            #                              dist_thres2,
            #                              rst2)
            s2 = None
        else:
            s2 = CountFromMask_Direct(ncell_from_inst)

        counter = Count_TwoStage(s1, s2)
    elif ty == 13:
        counter = Count_SumIntensity(50.)
    else:
        counter = Count_BlobDog(min_sigma=1., max_sigma=50., threshold=.5, exclude_border=False)
    return counter


class CountFromIntensityImage(ABC):
    @abstractmethod
    def count(self, im_cube: np.array) -> float:
        """
        Return the number of cells from image; The number should be counted one for each cell
        within the cube and 0.5 for each cell neighboring the boundaries (The counting of cells
        neighboring two or more walls of the cube area is 1/(num_neighbor + 1))

        The image is an float32 intensity image within range 0-1

        The return is a float estimate of (average of) the number of cells in the image given the image
        """
        raise NotImplementedError()

    @abstractmethod
    def interpretable(self, im_cube: np.array) -> np.array:
        """
        This method should return an overlay of the same size of the original image, intended for
        interpretation of the prediction result.
        """
        raise NotImplementedError()


class Count_SumIntensity(CountFromIntensityImage):
    def __init__(self, intensity_per_cell, thres=0.):
        """
        thres (float) - intensity below this is ignored; intensity above this is directly summed without change
        """
        self.intensity_per_cell = intensity_per_cell
        self.thres = thres
        assert intensity_per_cell > 0, f'{intensity_per_cell}'
        assert 0. <= thres <= 1., f'{thres}'

    def count(self, im_cube: np.array):
        mask = im_cube > self.thres
        ncells = (im_cube * mask).sum() / self.intensity_per_cell
        return ncells


class Count_BlobDog(CountFromIntensityImage):
    def __init__(self, **args):
        self.args = args

    def count(self, im_cube: np.array) -> float:
        blobs_dog = skimage.feature.blob_dog(np.array(im_cube, dtype=np.float32), **self.args)
        return blobs_dog.shape[0]

    def detect_blobs_draw(self, im_cube: np.array) -> np.array:
        blobs_dog = skimage.feature.blob_dog(np.array(im_cube, dtype=np.float32), **self.args)
        blobs_dog[:, 3:] = blobs_dog[:, 3:] * np.sqrt(3)  # adjust each sigma to get radius (this is still in pixels)

        blobs_img = np.zeros(im_cube.shape, dtype=np.float32)
        for i in range(blobs_dog.shape[0]):
            o = blobs_dog[i]
            iz, iy, ix, ir = int(o[0]), int(o[1]), int(o[2]), o[3]
            rr, cc = draw.circle_perimeter(iy, ix, int(max(1., ir / 2)))
            ylim, xlim = blobs_img.shape[1:3]
            inds = (rr >= ylim) + (cc >= xlim) + (rr < 0) + (cc < 0) == 0
            rr, cc = rr[inds], cc[inds]
            blobs_img[iz, rr, cc] = 1

        return blobs_img


class Count_TwoStage(CountFromIntensityImage):
    """
    First create a binary mask from input image, then use the input image and binary mask to
    get a cell num estimate
    """
    def __init__(self, s1: "MaskFromIntensityImage", s2: "CountFromMask"):
        self.s1 = s1
        self.s2 = s2

    def count(self, im_cube: np.array) -> float:
        mask = self.s1.mask(im_cube)
        ncells = self.s2.count(mask)
        return ncells


"""-------------------------Part 1: Mask from input image---------------------------------"""


class MaskFromIntensityImage(ABC):
    @abstractmethod
    def mask(self, im: np.array):
        """
        Create a binary mask from the input image
        """
        raise NotImplementedError()


class Mask_Thres(MaskFromIntensityImage):
    def __init__(self, thres):
        self.thres = thres

    def mask(self, im: np.array):
        return im > self.thres


class Mask_CellSeg3D(MaskFromIntensityImage):
    def __init__(self, cellseg3d_output):
        self.cellseg3d_output = cellseg3d_output

    def mask(self, im: np.array):
        return self.cellseg3d_output


"""-------------------------Part 2: cell count from mask---------------------------------"""


class CountFromMask(ABC):
    def __init__(self, ncell_from_inst: "NCellFromInst"):
        self.ncell_from_inst = ncell_from_inst

    def count(self, mask: np.array) -> float:
        lbl_im = self.inst(mask)
        return self.ncell_from_inst.count(lbl_im)

    def inst(self, mask: np.array) -> np.array:
        raise NotImplementedError()


class CountFromMask_Direct(CountFromMask):
    def inst(self, mask: np.array) -> np.array:
        lbl_im, nlbl = ndimage.label(mask)
        return lbl_im


class CountFromMask_Watershed(CountFromMask):
    def __init__(self, ncell_from_inst: "NCellFromInst", **args):
        super().__init__(ncell_from_inst)
        self.args = args

    def inst(self, mask: np.array) -> np.array:
        lbl_im = round_object_detection(mask, **self.args)
        return lbl_im


"""-------------------------Part 3: ncells from instance segmentation---------------------"""


class NCellFromInst(ABC):
    def count(self, inst: np.array) -> float:
        raise NotImplementedError()


class NCellFromInst_BySize(NCellFromInst):
    """
    Estimate the number of cells in cell clumps by looking at the size of the mask over it;
    small cells are still treated as one cell, while big clumps are counted as m cells where m
    depends on the number of voxels in the clump
    """

    def __init__(self, size_thres, voxel_per_cell):
        self.size_thres = size_thres
        self.voxel_per_cell = voxel_per_cell

    def count(self, inst: np.array) -> float:
        D, H, W = inst.shape
        nlbl = inst.max().item()
        mins = np.array([0, 0, 0], dtype=np.float32)
        maxs = np.array([D - 1, H - 1, W - 1], dtype=np.float32)
        ncells = 0.
        for i in range(1, nlbl + 1):
            inds = np.argwhere(inst == i)
            # get a 3-vector of number of bordering pixels in z, y, x dimensions
            nborders = np.sum((inds == mins[None, :]) + (inds == maxs[None, :]), axis=0)
            nborders = np.clip(nborders, 0, 1)
            nadd = 1 / (nborders.sum() + 1)
            if inds.shape[0] > self.size_thres:
                nadd += (inds.shape[0] - self.size_thres) / self.voxel_per_cell
            ncells += nadd
        return ncells


"""-------------------------Part 4: Implementations----------------------------------------"""


def instance_to_binary(seg_inst):
    """
    A convenient function to convert a [0, N] mask to binary mask while still somewhat preserving instance
    separation: Neighboring instances will have a thin dividing line drawn between them.
    """
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
    # `sure_bg = morph.binary_dilation(seg_bin, fp)`
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


"""-------------------------Part 5: Statistical Analysis----------------------------------"""


def Stats_PredictImages(counter: CountFromIntensityImage, im_list) -> list[float]:
    """
    Params
        counter (CountFromIntensityImage) - the counting algorithm used to count cells in images
        im_list (Iterable[float]) - a list of input images to be counted by counter algorithm
    Returns
        a list of float of estimated number of cells in each image
    """
    counted = []
    for im in im_list:
        counted_item = counter.count(im)
        counted.append(counted_item)
    return counted


def Stats_MAE(counted, gt):
    """
    Params
        counted (list) - the counted cells in each image
        gt (Iterable[float]) - the ground truth number of cells in each image
    Returns
        the mean absolute difference between counted and gt
    """
    return np.abs(np.array(counted, dtype=np.float32) - np.array(gt, dtype=np.float32)).mean().item()


def Stats_ShowScatterPairComparisons(counted: np.array, gt) -> None:
    """
    counted (a NCounter * NImages np.array) - the counted cells in each image, by each counter
    gt (Iterable[float]) - the ground truth number of cells in each image
    """
    gt_arr = np.array(gt, dtype=np.float32)
    fig, axes = plt.subplots(3, 6, figsize=(12, 6), sharex=True, sharey=True)
    for i in range(counted.shape[0]):
        X, Y = gt_arr, counted[i]

        iax, jax = i // 6, i % 6
        ax = axes[iax, jax]
        ax.set_box_aspect(1)
        ax.set_title(counting_method_inverse_dict[i])
        ax.scatter(X, Y)

