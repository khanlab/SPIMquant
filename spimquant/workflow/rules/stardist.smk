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
#    conda: '../envs/stardist.yaml'
    script: '../scripts/test_stardist.py'


