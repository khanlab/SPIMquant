import pandas as pd

low_abeta_subjects = (
    pd.read_csv("participants.tsv", sep="\t")
    .query(
        "group_label == 2 or group_label == 1 or (group_label ==3 and sex == 'Male')"
    )
    .participant_label.to_list()
)
subjects_by_group = dict()
spim_by_group=dict()

for group in ['1', '2', '3', '4']:
    subjects_by_group[group] = (
        pd.read_csv("participants.tsv", sep="\t")
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
            level=config["segment"]["fieldfrac_ds_level"],
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
            level=config["segment"]["fieldfrac_ds_level"],
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
            level=config["segment"]["fieldfrac_ds_level"],
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
     

rule deform_negative_mask_to_subject_nii:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        mask="negative_mask.nii",
        xfm_ras=rules.init_affine_reg.output.xfm_ras,
        invwarp=rules.deform_reg.output.invwarp,
    output:
        mask=bids(
            root=root,
            datatype="micr",
            desc="negative",
            level="{level}",
            from_="{template}",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN "
        "  -rm {input.mask} {output.mask} "
        "  -r {input.xfm_ras},-1 {input.invwarp}"
        #note: LABEL interpolation not possible with >1000 labels


