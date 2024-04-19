rule brainmask_penalty:
    """ generates a distance-based term to penalize blob detection
        at the brain boundary and outside it. Uses a signed distance transform
        (sdt) followed by scaled logistic function. 
    """ 
    input:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain',
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),  
    params:
        k=50, #steepness of logistic function (how fast it drops off) (penalty_weight)
        x0=0.1, #distance from boundary where it is penalized by 50%, in millimeters (penalty_distance_mm)
        

    output:
        sdt=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain',
            suffix="sdt.nii",
            **inputs["spim"].wildcards
        ),  
        penalty=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain',
            suffix="penalty.nii",
            **inputs["spim"].wildcards
        ),  
    shell:
        'c3d {input.mask} -sdt -scale -1 -o {output.sdt} -shift -{params.x0} -scale -{params.k} -exp -shift 1 -reciprocal '
        ' -o {output.penalty}'

rule blob_detection_betaamyloid:
    input:
        zarr=inputs["spim"].path,
        penalty=bids(
            root=root,
            datatype="micr",
            stain=config['masking']['stain'],
            level=config['masking']['level'],
            desc='brain',
            suffix="penalty.nii",
            **inputs["spim"].wildcards
        ),
    params:
        level=5,  #downsample-level to perform blob detection on
        min_sigma_um=1,
        max_sigma_um=100, # also serves as size of chunk borders
        threshold=0.06,
#        chunks=(1,400,100,100),
    output:
        npy=bids(
            root=root,
            datatype="micr",
            stain="{stain,BetaAmyloid}",
            suffix="blobs.npy",
            **inputs["spim"].wildcards
        ),
    script:
        '../scripts/blob_detection.py'


rule blob_detection_PI:
    input:
        zarr=inputs["spim"].path,
    params:
        level=3,  #downsample-level to perform blob detection on
        min_sigma_um=1,
        max_sigma_um=50, # also serves as size of chunk borders
        threshold=0.06,
        chunks=(1,400,200,200),
    output:
        npy=bids(
            root=root,
            datatype="micr",
            stain="{stain,PI}",
            suffix="blobs.npy",
            **inputs["spim"].wildcards
        ),
    script:
        '../scripts/blob_detection.py'


rule cellpose_BetaAmyloid:
    input:
        zarr=inputs["spim"].path,
    params:
        level=5,  #downsample-level to perform segmentation on
        chunks=(1,100,50,50),
    output:
        zarr=bids(
            root=root,
            datatype="micr",
            stain="{stain,BetaAmyloid}",
            desc='cellpose',
            suffix="dseg.zarr",
            **inputs["spim"].wildcards
        ),
    container: None
    script:
        '../scripts/cellpose.py'


"""
As a first pass, re-implement basic steps from following, but with skimage+dask:
  https://github.com/ChristophKirst/ClearMap/blob/master/ClearMap/ImageProcessing/SpotDetection.py
    
    Effectively this function performs the following steps:
        * illumination correction via :func:`~ClearMap.ImageProcessing.IlluminationCorrection.correctIllumination`
        * background removal via :func:`~ClearMap.ImageProcessing.BackgroundRemoval.removeBackground`
        * difference of Gaussians (DoG) filter via :func:`~ClearMap.ImageProcessing.Filter.filterDoG`
        * maxima detection via :func:`~ClearMap.ImageProcessing.MaximaDetection.findExtendedMaxima`
        * cell shape detection via :func:`~ClearMap.ImageProcessing.CellSizeDetection.detectCellShape`
        * cell intensity and size measurements via: :func:`~ClearMap.ImageProcessing.CellSizeDetection.findCellIntensity`,
          :func:`~ClearMap.ImageProcessing.CellSizeDetection.findCellSize`. 


 Parameters they used:


##############################################################################
# Image Processing
##############################################################################

correctIlluminationParameter = {
    "flatfield" : True,  # (str, True or None)  flat field intensities, if None d onot correct image for illumination, if True the 
    "background" : None, # (str, None or array) background image as file name or array, if None background is assumed to be zero
    "scaling" :  "Mean", # (str or None)        scale the corrected result by this factor, if 'max'/'mean' scale to keep max/mean invariant
    "save" : None,       # (str or None)        save the corrected image to file
    "verbose" : False    # (bool or int)        print / plot information about this step 
}

removeBackgroundParameter = {
    "size" : (15,15),  # size for the structure element of the morphological opening
    "save" : None,     # file name to save result of this operation
    "verbose" : False  # print / plot information about this step       
}


filterDoGParameter = {
    "size" :  (7, 7, 11),  # (tuple or None)      size for the DoG filter if None, do not correct for any background
    "sigma" : None,        # (tuple or None)      std of outer Guassian, if None autmatically determined from size
    "sigma2": None,        # (tuple or None)      std of inner Guassian, if None autmatically determined from size
    "save"  : None,        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose" : False      # (bool or int)        print / plot information about this step
}

findExtendedMaximaParameter = {
    "hMax" : 20,            # (float or None)     h parameter for the initial h-Max transform, if None, do not perform a h-max transform
    "size" : 5,             # (tuple)             size for the structure element for the local maxima filter
    "threshold" : 0,        # (float or None)     include only maxima larger than a threshold, if None keep all localmaxima
    "save"  : None,         # (str or None)       file name to save result of this operation if None dont save to file 
    "verbose" : False       # (bool or int)       print / plot information about this step
}

findIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "size" :  (3,3,3)       # (tuple)             size of the box on which to perform the *method*
}

detectCellShapeParameter = {
    "threshold" : 700,     # (float or None)      threshold to determine mask, pixel below this are background if None no mask is generated
    "save"  : None,        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose" : False      # (bool or int)        print / plot information about this step if None take intensities at the given pixels
}

"""


#perhaps illumination correction not needed if this is already done at the tile level

#background removal is via subtraction of morphological opening with disk (size=
#  can use: skimage.morphology.opening

#difference of gaussian
#  sigma=
#   can use skimage.feature.blob_dog


#maximum detection:
# can use: skimage.feature.peak_local_max
#

#cell shape detection:
#     imgws = watershed(-img, imgpeaks, mask = imgmask);
#
#     can use skimage.segmentation.watershed
#


#cell intensity and size measurements
# skimage.measure.regionprops
