import pandas as pd

low_abeta_subjects = (
    pd.read_csv("participants.tsv", sep="\t")
    .query(
        "group_label == 2 or group_label == 1 or (group_label ==3 and sex == 'Male')"
    )
    .participant_label.to_list()
)
print(1)
subjects_by_group = dict()
spim_by_group=dict()

for group in [1, 2, 3, 4]:
    print(2)
    subjects_by_group[group] = (
        pd.read_csv("participants.tsv", sep="\t")
        .query(f"group_label == {group}")
        .participant_label.to_list()
    )
    spim_by_group[group] = inputs["spim"].filter(subject=subjects_by_group[group])
    


print('Suspected Low Abeta Load subjects (used for generating neg mask)')
low_abeta_spim = inputs["spim"].filter(subject=low_abeta_subjects)
print(low_abeta_spim)

print(bids(root=root, stain="Abeta", group="lowload", suffix="avgfieldfrac.nii"))


print(subjects_by_group)




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
        lambda wildcards: spim_by_group[group].expand(
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
        avg_low_abeta_fieldfraction="avg_group-{group}_fieldfrac.nii",
    shell:
        "c3d {input} -mean -o {output}"
