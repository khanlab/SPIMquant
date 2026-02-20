# MRI registration

SPIMquant performs registration of MRI data, using BIDS to locate the MRI data for a subject, 
performing pre-processing on it (N4, motion-corrected averaging), registering with an MRI template (rigid+deformable)
to perform brain masking, then registration to the subject's SPIM dataset (affine+deformable).
Transformations between the MRI and the template (e.g. ABAv3) are then found by composing with warps from SPIM to template.

## Customizing snakebids filters to locate the MRI images:
By default, the snakebids filter to locate MRI images looks for a `T2w` suffix, as specified in the `pybids_inputs` in the `config/snakebids.yml`:
```
  mri:
    filters:
      suffix: 'T2w'
      extension: 'nii.gz'
      datatype: 'anat'
    wildcards:
      - subject
      - session
      - acquisition
      - run
      - reconstruction
      - suffix
      - extension
```
This can be overidden with the `--filter-mri` CLI option, e.g. with `--filter-mri suffix=T1w` to get T1w images. 

Depending on what is in your BIDS dataset, you may also need to specify additional filters, e.g. to point to the right session, acquisition  etc..

E.g. consider this more complicated example BIDS dataset snippet below:

```
sub-AS134F1
├── ses-11m
│   └── anat
│       ├── sub-AS134F1_ses-11m_acq-GRE_rec-avgecho_run-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_rec-avgecho_run-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_rec-avgecho_run-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_rec-avgecho_run-4_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_rec-avgecho_run-5_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_rec-avgecho_run-6_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-1_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-1_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-1_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-2_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-2_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-2_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-3_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-3_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-3_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-4_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-4_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-4_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-5_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-5_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-5_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-6_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-6_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-11m_acq-GRE_run-6_echo-3_T2starw.nii.gz
│       └── sub-AS134F1_ses-11m_acq-TOF_run-1_angio.nii.gz
├── ses-14m
│   └── anat
│       ├── sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-4_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-5_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-6_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-1_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-1_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-1_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-2_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-2_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-2_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-3_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-3_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-3_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-4_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-4_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-4_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-5_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-5_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-5_echo-3_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-6_echo-1_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-6_echo-2_T2starw.nii.gz
│       ├── sub-AS134F1_ses-14m_acq-GRE_run-6_echo-3_T2starw.nii.gz
│       └── sub-AS134F1_ses-14m_acq-TOF_run-1_angio.nii.gz
└── ses-exvivo14m1
    └── micr
        ├── sub-AS134F1_ses-exvivo14m1_sample-brain_acq-imaris4x_SPIM.json
        └── sub-AS134F1_ses-exvivo14m1_sample-brain_acq-imaris4x_SPIM.ome.zarr
```

Here we have multiple MRI sessions, and multiple T2starw anatomicals, multiple runs of 
individual echoes, and multiple runs of average echo (rec-avgecho).

In this case we would want to import the rec=avgecho runs, from ses=14m (closest in time
to the exvivo session where the microscopy lives). This is done with the following CLI call:

```
pixi run spimquant /nfs/trident3/mri/prado/ki3/bids output_dir participant \
   --register-to-mri --no-segmentation
   --filter-mri suffix=T2starw reconstruction=avgecho session=14m \
   -c all  --mri-resample-percent 200
```

This ends up selecting these images for the registration:
```
sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-1_T2starw.nii.gz
sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-2_T2starw.nii.gz
sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-3_T2starw.nii.gz
sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-4_T2starw.nii.gz
sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-5_T2starw.nii.gz
sub-AS134F1_ses-14m_acq-GRE_rec-avgecho_run-6_T2starw.nii.gz
```
which the workflow will use to produce a single motion-corrected average. 

Note, we are also running the workflow without SPIM segmentation `--no-segmentation` to skip 
the computationally-intensive SPIM quantification and just do registration. 
We are also using the `--mri-resample-percent` option to use an upsampled reference when performing
the motion-corrected averaging, which provides some super-resolution enhancement. 
