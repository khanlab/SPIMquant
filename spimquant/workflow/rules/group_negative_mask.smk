# THESE ARE DATASET-SPECIFIC RULES, TO BE REMOVED LATER..

import pandas as pd

low_abeta_subjects = (
    pd.read_csv("wip_mouse/participants.tsv", sep="\t")
    .query(
        "group_label == 2 or group_label == 1 or (group_label ==3 and sex == 'Male')"
    )
    .participant_label.to_list()
)
subjects_by_group = dict()
spim_by_group=dict()

for group in ['1', '2', '3', '4']:
    subjects_by_group[group] = (
        pd.read_csv("wip_mouse/participants.tsv", sep="\t")
        .query(f"group_label == {group}")
        .participant_label.to_list()
    )
    spim_by_group[group] = inputs["spim"].filter(subject=subjects_by_group[group])
    

print('first group')
print(spim_by_group['1'])
#print('Suspected Low Abeta Load subjects (used for generating neg mask)')
low_abeta_spim = inputs["spim"].filter(subject=low_abeta_subjects)
#print(low_abeta_spim)

#print(bids(root=root, stain="Abeta", group="lowload", suffix="avgfieldfrac.nii"))


#print(subjects_by_group)


#print(spim_by_group)

rule avg_fieldfrac_low_abeta:
    input:
        low_abeta_spim.expand(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{level}",
                desc="{desc}",
                space="{template}",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards
            ),
            stain="Abeta",
            level=config["downsampling_level"],
            desc="otsu",
            template=config["template"],
        ),
    output:
        avg_low_abeta_fieldfraction="avg_lowload_Abeta_fieldfrac.nii",
    shell:
        "c3d {input} -mean -o {output}"

rule avg_fieldfrac_bygroup:
    input:
         lambda wildcards: spim_by_group[wildcards.group].expand(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{level}",
                desc="{desc}",
                space="{template}",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards
            ),
            stain="Abeta",
            level=config["downsampling_level"],
            desc="otsu",
            template=config["template"],
        ),
    output:
        avg="avg_group-{group}_fieldfrac.nii",
    shell:
        "c3d {input} -mean -o {output}"

rule avg_masked_fieldfrac_bygroup:
    input:
         lambda wildcards: spim_by_group[wildcards.group].expand(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{level}",
                desc="{desc}",
                space="{template}",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards
            ),
            stain="Abeta",
            level=config["downsampling_level"],
            desc="otsupenalty",
            template=config["template"],
        ),
    output:
        avg="avg_group-{group}_maskedfieldfrac.nii",
    shell:
        "c3d {input} -mean -o {output}"



rule create_negative_mask:
    input:
        avg_low_abeta_fieldfraction="avg_lowload_Abeta_fieldfrac.nii",
    output:
        negative_mask="negative_mask.nii"
    shell:
        'c3d {input} -smooth 3x3x3vox  -threshold 2 inf 0 1 -o {output}'
     

rule avg_roi_fieldfrac_bygroup:
    input:
         lambda wildcards: spim_by_group[wildcards.group].expand(
        bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            space="{template}",
            stain="{stain}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards
        ),
            stain="Abeta",
            template=config["template"],
            seg=wildcards.seg,
        ),
    output:
        avg="avg_seg-{seg}_group-{group}_fieldfrac.nii",
    shell:
        "c3d {input} -mean -o {output}"


