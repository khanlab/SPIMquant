rule get_model:
    input:
        model_zip=lambda wildcards: storage(config['stardist']['model_urls'][wildcards.model_name])
    output:
        model_dir=directory('resources/stardist_models/{model_name}')
    shell:
        'unzip -d {output.model_dir} {input.model_zip}'
        
rule test_stardist:
    input:
        model_dir=directory('resources/stardist_models/3D_demo')
    output:
        out_dir=directory('test_stardist')
    script: '../scripts/test_stardist.py'

rule segment_stardist:
    input:
        spim=inputs["spim"].path,
    params:
        zarrnii_kwargs=lambda wildcards: {"orientation": config["orientation"],channel_labels=[wildcards.stain]},
    output:
        regionprops=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="correctedn4",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        biasfield=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="n4",
                    suffix="biasfield.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=500000,
        disk_mb=2097152,
        runtime=60,
    script:
        "../scripts/n4_biasfield.py"



