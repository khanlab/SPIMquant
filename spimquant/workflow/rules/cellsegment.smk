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


cellseg3d_dataset_name = config["cellsegment"]["init_dataset"]["dataset_name"]
zarr_path_expanded = cconfig.inputs['spim'].expand()[0]
rule cellseg3d_preprocess:
    """
    Training set of CellSeg3D requires to sample slices from the OME_zarr file.
    This command creates the required slices.
    """
    params:
        zarr = zarr_path_expanded,
        output_dir=config["output_dir"],
        cellsegment=config["cellsegment"]["init_dataset"],  # pass in the dictionary
        command='init_dataset',
    output:
        dataset_dir=directory(f'{config["output_dir"]}/{cellseg3d_dataset_name}')
    script:
        "../../../../spimquant_CellSeg3D/spimquant.py"

rule cellseg3d_train:
    """
    Train CellSeg3D model on the slices output by cellseg3d_train rule.
    This command creates a model config folder that contains the trained model.
    """
    input:
        dataset_dir=f'{config["output_dir"]}/{cellseg3d_dataset_name}'
    params:
        output_dir=config["output_dir"],
        cellsegment=config["cellsegment"]["init_dataset"] | config["cellsegment"]["train"],
        command='train',
    output:
        model_config=directory(f'{config["output_dir"]}/{cellseg3d_dataset_name}_model_config')
    container:
        None  # TODO: put this in container
    script:
        "../../../../spimquant_CellSeg3D/spimquant.py"

rule cellseg3d_predict:
    """
    Taking a model config folder, predicting on a new 3d numpy chunk of arbitrary size that is >= 64 * 64 * 64
    The output of this rule is a .npy file containing the normalized probabilities (0-1) for each classes among the 
    num_classes classes in the intermediate representation.
    """
    input:
        model_config=f'{config["output_dir"]}/{cellseg3d_dataset_name}_model_config'
    params:
        zarr = zarr_path_expanded,
        output_dir=config["output_dir"],
        cellsegment=config["cellsegment"]["init_dataset"] | config["cellsegment"]["train"] | config["cellsegment"]["predict"],
        command='predict',
    output:
        pred_result_dir=directory(f'{config["output_dir"]}/{cellseg3d_dataset_name}_pred')
    container:
        None  # TODO: put this in container
    script:
        "../../../../spimquant_CellSeg3D/spimquant.py"

rule cellseg3d_view:
    """
    View the prediction result with Napari; overlay the results on original image
    The result is a nclass*Depth*Height*Width float image
    """
    input:
        pred_result_dir=f'{config["output_dir"]}/{cellseg3d_dataset_name}_pred'
    script:
        '../scripts/cellseg3d_napari_pair_overlay.py'

rule cellseg3d_supervised_dataset_gen:
    """
    Generate supervised datasets for annotations
    """
    input:
        model_config=f'D:/progtools/RobartsResearch/data/lightsheet/mousebrain_chan0_20240627_1_model_config'
    params:
        zarr=zarr_path_expanded,
        cellsegment=config["cellsegment"]["init_dataset"] | config["cellsegment"]["train"] | config["cellsegment"][
            "predict"],
    output:
        dataset_dir=directory(f'{config["output_dir"]}/supervised_datasets')
    script:
        "../../../../spimquant_CellSeg3D/supervised_dataset_gen.py"
