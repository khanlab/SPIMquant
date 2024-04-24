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
            desc="brain",
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),
    params:
        k=50,  #steepness of logistic function (how fast it drops off) (penalty_weight)
        x0=0.1,  #distance from boundary where it is penalized by 50%, in millimeters (penalty_distance_mm)
    output:
        sdt=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="sdt.nii",
            **inputs["spim"].wildcards
        ),
        penalty=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="penalty.nii",
            **inputs["spim"].wildcards
        ),
    shell:
        "c3d {input.mask} -sdt -scale -1 -o {output.sdt} -shift -{params.x0} -scale -{params.k} -exp -shift 1 -reciprocal "
        " -o {output.penalty}"


rule blob_detection_betaamyloid:
    input:
        zarr=inputs["spim"].path,
        penalty=bids(
            root=root,
            datatype="micr",
            stain=config["masking"]["stain"],
            level=config["masking"]["level"],
            desc="brain",
            suffix="penalty.nii",
            **inputs["spim"].wildcards
        ),
    params:
        level=lambda wildcards: int(wildcards.level),  #downsample-level to perform blob detection on
        min_sigma_um=1,
        max_sigma_um=100,  # also serves as size of chunk borders
        threshold=0.06,
        chunks=(1,127,122,116),
    output:
        sparse_npz=bids(
            root=root,
            datatype="micr",
            level='{level}',
            stain="{stain,BetaAmyloid}",
            suffix="sparseblobs.npz",
            **inputs["spim"].wildcards
        ),
        points_npy=bids(
            root=root,
            datatype="micr",
            level='{level}',
            stain="{stain,BetaAmyloid}",
            suffix="points.npy",
            **inputs["spim"].wildcards
        ),
    threads: 6
    container: None
    shadow: 'minimal'
    script:
        "../scripts/blob_detection.py"

rule filter_blobs_betaamyloid:
    input:
        points_npy=bids(
            root=root,
            datatype="micr",
            level='{level}',
            stain="{stain,BetaAmyloid}",
            suffix="points.npy",
            **inputs["spim"].wildcards
        ),
        penalty=bids(
            root=root,
            datatype="micr",
            stain=config["masking"]["stain"],
            level=config["masking"]["level"],
            desc="brain",
            suffix="penalty.nii",
            **inputs["spim"].wildcards
        ), 
    params:
        #penalty is from 0 to 1 (0.5 at the penalty_distance)
        #  could consider simply using the distance map directly, instead of the sigmoid-filtered one.. 
        #  - the sigmoid penalty is better suited for when a soft-penalty is needed, not a hard threshold
        threshold=0.5
    output:
        points_npy=bids(
            root=root,
            datatype="micr",
            level='{level}',
            desc='filtered',
            stain="{stain,BetaAmyloid}",
            suffix="points.npy",
            **inputs["spim"].wildcards
        ),
    script:
        '../scripts/filter_blobs.py'

rule map_labels_to_blobs:
    input:
        points_npy=bids(
            root=root,
            datatype="micr",
            level=config['blobdetect']['level'],
            desc='filtered',
            stain="{stain}",
            suffix="points.npy",
            **inputs["spim"].wildcards
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            desc="deform",
            level=config['blobdetect']['dseg_level'],
            from_=config['blobdetect']['dseg_template'],
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
        ),
    output:
        cells_tsv=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            suffix="blobs.tsv",
            **inputs["spim"].wildcards
        ),
    script: 
        '../scripts/map_labels_to_blobs.py'


rule generate_subject_density_tsv:
    """ this reads in the dseg to calc volume, blobs.tsv to calculate number per label, and calculates the density, writing out density and volume for each label in a new tsv file (number of rows = number of labels in atlas)"""
    input:
        dseg=bids(
            root=root,
            datatype="micr",
            desc="deform",
            level=config['blobdetect']['dseg_level'],
            from_=config['blobdetect']['dseg_template'],
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
        ),
        cells_tsv=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            suffix="blobs.tsv",
            **inputs["spim"].wildcards
        ),
    output:
        density_tsv=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            suffix="blobdensity.tsv",
            **inputs["spim"].wildcards
        ),
    script: '../scripts/generate_subject_density_tsv.py'



rule generate_subject_density_nii:
    """ generates a heatmap """
    input:
         density_tsv=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            suffix="blobdensity.tsv",
            **inputs["spim"].wildcards
        ),   
        dseg=bids(
            root=root,
            datatype="micr",
            desc="deform",
            level=config['blobdetect']['dseg_level'],
            from_=config['blobdetect']['dseg_template'],
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
        ),
    output:
         density_nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            suffix="blobdensity.nii",
            **inputs["spim"].wildcards
        ),
    script: '../scripts/generate_subject_density_nii.py'
