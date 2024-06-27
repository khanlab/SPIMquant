rule cellpose:
    """ haven't tried on PI data or full-res, but so far doesn't really work well.."""
    input:
        zarr=cconfig.inputs["spim"].path,
    params:
        level=lambda wildcards: int(wildcards.level),  #downsample-level to perform segmentation on
        chunks=(200, 200, 200),
    output:
        zarr=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="cellpose",
                suffix="dseg.zarr",
                **cconfig.inputs["spim"].wildcards
            )
        ),
    container:
        None  # TODO: not in container yet, put this in container
    threads: 6
    script:
        "../scripts/cellpose.py"


checkpoint cellseg3d_create_trainset:
    """
    Training set of CellSeg3D requires to sample slices from the OME_zarr file.
    This command creates the required slices.
    """
    params:
        zarr = 'gs://khanlab-lightsheet/data/onuska_lifecanvas/bids/sub-o21/micr/sub-o21_sample-brain_acq-prestitched_SPIM.ome.zarr',
        output_dir=config["output_dir"],
        cellsegment=config["cellsegment"],  # pass in the dictionary
        command='init_dataset',
    output:
        dataset_dir=directory(f'{config["output_dir"]}/{config["cellsegment"]["dataset_name"]}')
    script:
        "../../../../spimquant_CellSeg3D/spimquant.py"

rule cellseg3d_train:
    """
    Train CellSeg3D model on the slices output by cellseg3d_train rule.
    This command creates a model config folder that contains the trained model.
    """
    input:
        dataset_dir=f'{config["output_dir"]}/{config["cellsegment"]["dataset_name"]}'
    params:
        output_dir=config["output_dir"],
        cellsegment=config["cellsegment"],# pass in the dictionary
        command='train',
    output:
        model_config=directory(f'{config["output_dir"]}/{config["cellsegment"]["dataset_name"]}_model_config')
    container:
        None  # TODO: not in container yet, put this in container
    script:
        "../../../../spimquant_CellSeg3D/spimquant.py"

rule cellseg3d_predict:
    """
    Taking a model config folder, predicting on a new numpy chunk of arbitrary size > 64*64*64
    """
    input:
        model_config=f'{config["output_dir"]}/{config["cellsegment"]["dataset_name"]}_model_config'
    params:
        zarr = 'gs://khanlab-lightsheet/data/onuska_lifecanvas/bids/sub-o21/micr/sub-o21_sample-brain_acq-prestitched_SPIM.ome.zarr',
        output_dir=config["output_dir"],
        cellsegment=config["cellsegment"],# pass in the dictionary
        command='predict',
    output:
        pred_result_dir=directory(f'{config["output_dir"]}/{config["cellsegment"]["dataset_name"]}_pred')
    container:
        None  # TODO: not in container yet, put this in container
    script:
        "../../../../spimquant_CellSeg3D/spimquant.py"
