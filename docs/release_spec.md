# SPIMquant release contract

**Schema version:** 0.1.0 (draft) · **Status:** proposed, not yet implemented

A *release* is the exact set of files SPIMquant hands to MIND (the data-management
portal and the public atlas-space explorer). Nothing outside the release is uploaded,
and nothing in the release is a work file. This page fixes the interface that two
pieces of work share:

- **PR1, the `release` target:** writes the participant and group releases described
  here, and a CI check validates them against this schema on the test dataset.
- **PR2, `--group-manifest`:** group-level analysis reads subjects and covariates from
  an aggregation manifest, and stages `mind://` subjects by fetching release tables only.

MRI pipelines are not part of SPIMquant, but their derivatives follow the
[MRI derivative schema](#mri-derivative-schema) below so the portal and explorer can
treat them like SPIM releases.

Each section is built on what SPIMquant writes today (see
[Output Files Reference](reference/outputs.md)). Differences from current outputs are
listed in [Changes needed in SPIMquant](#changes-needed-in-spimquant).

---

## Release levels

| Level | Made by | When | Contents |
| --- | --- | --- | --- |
| **Participant release** | `release` target at the participant level, one subject (and session) per run | On ingestion into MIND, or by any participant-level run | Everything about one subject: tables, view copies, sidecars, provenance. Nothing that needs other subjects. |
| **Group release** | `release` target at the group level, from an aggregation manifest | Later, over participant releases selected across datasets or projects | Aggregated tables, group maps, covariates, the release manifest (explorer entry point, DOI). |

**Participant releases are the unit.** MIND runs SPIMquant on each subject as it is
ingested, so most participant releases are made one subject at a time and never have
group-level outputs or metadata of their own. Group releases (and group analysis) are
built afterwards by aggregating participant releases with an
[aggregation manifest](#aggregation-manifest-group-manifest), usually the ingestion
releases stored in MIND. So:

- A participant release is complete on its own: it depends on no other subject, no
  dataset-level table (e.g. `participants.tsv`) and no group-level output.
- **Covariates are not part of a participant release.** They come from the dataset's
  `participants.tsv` (or the MIND catalog) when an aggregation manifest is built.
- **QC status is not part of a participant release** either: it is recorded in the
  MIND catalog after review, and read when a group release is built.
- **Each run is its own derivative dataset in MIND.** A run's release (one subject on
  ingestion, or several) is stored as a separate derivative dataset, with its own
  `dataset_description.json` and label tables; releases from different runs are never
  merged into one dataset. Aggregation manifests point at the subject directories
  inside these datasets, and the compatibility check compares their provenance.
- Whether a participant release was made alone or in a multi-subject run makes no
  difference to its subject directory.

---

## Conventions

| Item | Rule |
| --- | --- |
| Template | `ABAv3` for every release that MIND accepts. Template-space data and template coordinates are in ABAv3 physical space. |
| Label table | `tpl-ABAv3_seg-<seg>_dseg.tsv` from the template's Zenodo record (columns `index`, `name`, ...). The release copies the label tables it uses, so a release is self-contained. |
| Canonical labels | `seg-all` (the full ABAv3 label set). Instance tables carry `seg-all` labels only. Coarser parcellations come from the label hierarchy or from the regional table. |
| Released parcellations | Regional tables cover the template's `default_segs` (`all`, `coarse`, `mid`, `fine` for ABAv3). |
| Units | Coordinates and distances in **mm**; volumes in **mm³**; densities in **objects/mm³**; field fraction in **percent (0–100)**; intensities in arbitrary units. Units are written in the table (`units` column) or in the JSON sidecar, never assumed. |
| Subject key | `participant_id` (`sub-<label>`) and `session_id` (`ses-<label>`, or null when the dataset has no sessions). Every table carries both, so tables from many subjects concatenate without joins. |
| Scans | A participant release holds every scan of a subject and session: one **reference** scan and any **additional** scans (`--aux_spim`, e.g. a zoomed-in field of view), told apart by `acq-`. Every table carries `acq` and `scan_role` (`reference` or `additional`); file names carry `acq-` as usual. See [Scans](#scans). |
| Other BIDS entities | `sample` and similar appear in file names and in the provenance sidecar. |
| Coordinates | Every image and coordinate column states its coordinate system; view copies, NIfTI files and table coordinates of the same space agree. See [Coordinate systems](#coordinate-systems). |
| Sample and acquisition metadata | Antibodies, fluorophores, channel optics, microscope and sample preparation, in each scan's BIDS sidecar. See [Sample and acquisition metadata](#sample-and-acquisition-metadata). |
| Stain / channel | The stain name as in the BIDS `stain-` entity (e.g. `Abeta`, `Iba1`, `CD31`). Tables use the generic column name `channel` so SPIM and MRI rows share a schema; for SPIM rows `channel` is the stain. |
| Segmentation method | Two columns in every table: `seg_method` (mask method, e.g. `lantern`, `otsu+k3i2`) and `object_method` (`components`, `watershed`). Per-stain files such as masks carry them in the file name (`desc-`, `obj-`); results and release files carry neither. See [Segmentation methods per stain](#segmentation-methods-per-stain). |
| Tables | Parquet (zstd compression) for everything the explorer reads. TSV is kept only where SPIMquant already writes it for users. Every table has a JSON sidecar describing its columns (see `create_json_sidecar.py`). |
| Missing values | Null, never `n/a` or `0`. A region with no detected objects has `count = 0`; a region outside the field of view has null metrics. |
| Versioning | `schema_version` is semver. A major bump means a release reader for the old version cannot read the new one (renamed or removed column, changed units). Adding optional columns or metrics is a minor bump. |

### Typical MIND stain panel

The design case is **Abeta + Iba1 + CD31**:

| Stain | What SPIMquant does | Where it lands in the release |
| --- | --- | --- |
| Abeta | Plaques: `lantern` masks, split into objects by `components` (built-in default) | Regional table rows; instance table rows (one per plaque); `Abeta_dist` and `Abeta_edgedist` on every object |
| Iba1 | Microglia: `gmm+n3k1` masks by default; `--object_method Iba1=watershed` so touching cells count separately | Regional table rows; instance table rows (one per cell); `Iba1_dist` and `Iba1_edgedist` on every object |
| CD31 | Vessels: `vesselfm` mask and its signed distance transform (SDT) | `CD31_sdt` on every Abeta and Iba1 object, and its regional mean; template-space vessel field-fraction map. The vessel graph (with calibre per edge) becomes the CD31 instance table later. |

Abeta and Iba1 objects are colocalized with each other even though they use different
methods.

---

## Scans

A sample can be imaged more than once, e.g. a whole brain at 1x and a 12x field of view
of a region. SPIMquant registers the **reference** scan to the template, and each
**additional** scan (`--aux_spim`) to the reference scan; the additional scan gets atlas
labels and template coordinates through it.

- Both kinds of scan are in the same participant release, with their own tables, view
  copies and sidecars, named with their `acq-` entity.
- Tables carry `acq` (null when the dataset has no `acq-` entity) and `scan_role`.
- The regional table of an additional scan has the `coverage` metric: the fraction of
  each region (as labelled in the reference scan) inside its field of view. Volumes,
  counts and densities of an additional scan are over the covered part only.
- A subject's scans are **never pooled**: group tables, group statistics and the
  explorer keep each acquisition separate, as SPIMquant's group level does today, so a
  subject never counts twice. Group statistics of additional scans drop regions they
  cover less than `--aux_min_coverage` of.
- An additional scan's view copy comes with its transform to the reference scan, so a
  viewer can place the zoomed field of view inside the reference volume (see
  [Coordinate systems](#coordinate-systems)).

---

## Coordinate systems

Overlays only line up if every image and every coordinate column of a space use the
same frame. OME-Zarr 0.5 transforms hold only a per-axis scale and translation (no
rotation or flip; general transforms are still an open OME-Zarr proposal), while NIfTI
affines can hold any affine. So the release fixes these rules:

| Space | Frame | Images | Coordinate columns |
| --- | --- | --- | --- |
| Template (`space-ABAv3`) | ABAv3 world coordinates in mm, as defined by the template's NIfTI affine | NIfTI as SPIMquant writes them; OME-Zarr view copies on the same grid, axes `z`, `y`, `x`, with `scale` and `translation` reproducing the NIfTI affine | `template_x`, `template_y`, `template_z` |
| Subject (each scan) | The scan's physical frame in mm, as zarrnii reads it (with the input's orientation applied) | OME-Zarr view copies with the same `scale` and `translation`; anatomical orientation in the `xyz_orientation` attribute | `pos_x`, `pos_y`, `pos_z` |

- **Axis-aligned only.** A frame that needs a flip or rotation to match is written with
  the flip applied to the voxel order, so `scale` stays positive and the transform stays
  scale-and-translation. If the ABAv3 NIfTI affine is not diagonal, template-space view
  copies are resampled to a diagonal grid of the same world space.
- **Axis order and units.** OME-Zarr axes are `z`, `y`, `x` (with `c` for channels when
  present), unit `millimeter`. Table coordinate columns are named for the world axes
  (`x`, `y`, `z`), whatever the array order.
- **Additional scans.** The transform from an additional scan to its reference scan
  (from `--aux_spim`) is released as `xfm/<ent>_acq-<acq>_from-<acq>_to-<ref acq>_xfm.json`:
  a 4×4 affine in mm, world-to-world, in the axis order `x`, `y`, `z`, which a viewer
  applies as a layer transform. The ITK/greedy text file SPIMquant writes is released
  alongside it.
- **Validation.** For a sample of objects, the voxel their table coordinates fall in, in
  every view copy of that space, must hold the expected label or mask value (see
  [Validation](#validation)).

---

## Sample and acquisition metadata

What was stained with what, and how it was imaged, lives in each scan's BIDS
microscopy sidecar, `micr/<ent>_SPIM.json`. SPIMquant does not create this metadata:
the release copies the input sidecar unchanged, as SPIMpack wrote it at import.

**BIDS microscopy fields** (used as defined by BIDS):

| Group | Fields |
| --- | --- |
| Microscope | `Manufacturer`, `ManufacturersModelName`, `DeviceSerialNumber`, `StationName`, `SoftwareVersions` |
| Institution | `InstitutionName`, `InstitutionAddress`, `InstitutionalDepartmentName` |
| Acquisition | `PixelSize`, `PixelSizeUnits` (required), `Immersion`, `NumericalAperture`, `Magnification`, `ImageAcquisitionProtocol`, `OtherAcquisitionParameters` |
| Sample | `BodyPart`, `BodyPartDetails`, `SampleEnvironment`, `SampleFixation`, `SampleEmbedding`, `SliceThickness`, `TissueDeformationScaling`, `SampleExtractionProtocol` |
| Per channel | `SampleStaining`, `SamplePrimaryAntibody`, `SampleSecondaryAntibody` (antibodies as an RRID, or name, supplier and catalog number) |

**Extension fields** for per-channel optics, which BIDS has no fields for
(`OtherAcquisitionParameters` is a single free-text string). They follow the BIDS
convention for per-channel values:

| Field | Units | Meaning |
| --- | --- | --- |
| `SampleFluorophore` | | Fluorophore or dye of each channel |
| `ExcitationWavelength` | nm | Excitation wavelength of each channel |
| `EmissionWavelength` | nm | Centre of the emission filter of each channel |
| `LaserPower` | percent | Laser power setting of each channel |
| `ExposureTime` | ms | Camera exposure of each channel |

**Per-channel rule.** A per-channel field is either a single value (the same for every
channel) or a list with one value per channel, in the order of `SampleStaining`, which
is the order of the image's channels. Release validation fails if a list has a
different length from `SampleStaining`, or if `SampleStaining` does not match the
scan's stains.

The release manifest summarises these per stain (antibodies, fluorophore) and per scan
(microscope, magnification); see
[Release manifest](#release-manifest-release_manifestjson). What a stain marks
(`Target`, e.g. `amyloid plaques`) is given when the release is built.

SPIMpack writes the BIDS fields at import; the extension fields need adding to its
sidecar field list. The `.ims` metadata holds most of the channel optics, so SPIMpack can
fill them in without asking the user.

---

## Segmentation methods per stain

Results combine stains by **stain**, not by method: each segmented stain has one
*primary* mask method and one primary object method, and every cross-stain result
(colocalization, merged regional table, instance table, distances, group tables) uses
the primaries. Stains can use different methods and still be compared.

```
--seg_method SPEC [SPEC ...]
    STAIN=METHOD[,METHOD...]   per stain; the first is primary, the rest are comparisons
    METHOD                      every stain not named
    METHOD: th<v> | otsu+k<k>i<i> | gmm+n<n>k<k> | dog+s<s>r<r>k<k> | lantern
    built-in defaults: Abeta=lantern (also abeta, BetaAmyloid) when a GPU is available,
                       else gmm+n3k1; every other stain gmm+n3k1

--object_method SPEC [SPEC ...]
    STAIN=OBJ[,OBJ...] | OBJ   same rules as --seg_method
    OBJ: watershed | components | none
    built-in default: components for every stain

--gpu auto|yes|no                 auto: torch built with CUDA (the gpu pixi env)

--stains_for_seg STAIN [...]   which stains are segmented; stains named above are added
```

- **Precedence:** `STAIN=` on the command line, then a bare value on the command line,
  then the built-in per-stain default, then the built-in global default. The resolved
  plan is printed at start-up and recorded in the provenance sidecar.
- **Object methods:** `components` is connected components (touching objects merge);
  `watershed` is a seeded watershed of the mask, seeded at peaks of the smoothed
  bias-corrected image (today's `+inst`); `none` skips every object-level output for
  that stain (counts, density, instance rows, colocalization, distances), leaving field
  fraction and volume.
- **GPU methods** (`lantern`, future models) have a CPU fallback in the config
  (`gpu_seg_method_fallbacks`). Without a GPU a default GPU method falls back,
  with a note; one given on the command line is an error.
- **`+inst` is rejected** with an error pointing to `--object_method watershed`.
- **Comparison methods** (a second method for a stain) produce per-stain outputs and
  extra rows in the regional table. They are not in the instance table,
  colocalization or group statistics.
- **Distances are automatic:**
  - for every stain with objects, each object gets `<stain>_dist` (centroid to nearest
    centroid of that stain, mm; for its own stain the nearest *other* object) and
    `<stain>_edgedist` (`<stain>_dist` minus the nearest object's equivalent-sphere
    radius from its volume, clipped at 0; an approximate distance to its surface);
  - when vessels are run, each object gets `<vessel stain>_sdt` (signed distance to the
    vessel mask, mm, negative inside), sampled from the SDT the vessel graph already
    needs. SDTs of other stains are off; a config key (`sdt_stains`) can enable them.
- **File names** (snakebids orders the entities `stain`, `space`, `seg`, `from`,
  `level`, `obj`, `desc`): per-stain files carry the methods, `desc-<seg_method>` on
  masks and mask-derived maps (`<ent>_stain-Abeta_level-0_desc-lantern_mask.ozx`) and
  `obj-<object_method>` before it on per-stain object outputs
  (`<ent>_stain-Iba1_seg-all_from-ABAv3_level-5_obj-watershed_desc-gmm+n3k1_segstats.tsv`).
  Cross-stain results and release files carry neither
  (`<ent>_seg-all_from-ABAv3_mergedsegstats.tsv`, `<ent>_space-ABAv3_regionprops.parquet`);
  the methods are in table columns and sidecars.
- Vessels keep their own options (`--stains_for_vessels`, `--no_vessels`) and outputs
  (mask, SDT, graph); they have no object method.

---

## Participant release

### Layout

The release is a BIDS-derivative dataset of its own, written under the release root
(`<output_dir>/release/` by default, or a `mind://` URI; see
[Release root and storage format](#release-root-and-storage-format)).
A run of one subject writes one subject directory; a run of many writes one each, with
the same files.
Paths below are relative to that root; `<ent>` is a scan's entity prefix
(`sub-01_ses-01_sample-brain_acq-1x`, without `stain-`). Each scan of the subject
(reference and additional) has the per-scan files below.

```
release/
├── dataset_description.json              # DatasetType: derivative, GeneratedBy: spimquant (this run)
├── tpl-ABAv3/
│   ├── tpl-ABAv3_seg-<seg>_dseg.tsv      # label tables used by this release
│   └── tpl-ABAv3_seg-<seg>_dseg.json
└── sub-<label>/[ses-<label>/]
    ├── sub-<label>[_ses-<label>]_provenance.json   # provenance sidecar, all scans (below)
    ├── micr/
    │   └── <ent>_SPIM.json                             # input BIDS microscopy sidecar, copied
    ├── xfm/                                            # additional scans only
    │   ├── <ent>_from-<acq>_to-<ref acq>_xfm.json      # 4x4 affine to the reference scan
    │   └── <ent>_from-<acq>_to-<ref acq>_xfm.txt       # same, as written by SPIMquant
    ├── tabular/
    │   ├── <ent>_space-ABAv3_roistats.parquet          # regional table, long (all segs, stains, descs)
    │   ├── <ent>_space-ABAv3_roistats.json
    │   ├── <ent>_seg-<seg>_from-ABAv3_mergedsegstats.tsv   # wide TSV, primary methods, for users
    │   ├── <ent>_seg-<seg>_from-ABAv3_mergedsegstats.json
    │   ├── <ent>_space-ABAv3_instances.parquet         # instance table (all stains, primary methods)
    │   ├── <ent>_space-ABAv3_instances.json
    │   ├── <ent>_space-ABAv3_coloc.parquet             # colocalization pairs (2+ stains with objects)
    │   └── <ent>_space-ABAv3_coloc.json
    ├── view/                                           # OME-Zarr view copies
    │   ├── <ent>_stain-<stain>_SPIM.ome.zarr            # subject space, full resolution pyramid
    │   ├── <ent>_stain-<stain>_space-ABAv3_SPIM.ome.zarr
    │   ├── <ent>_stain-<stain>_mask.ome.zarr            # subject-space primary masks
    │   └── <ent>_seg-all_from-ABAv3_dseg.ome.zarr       # subject-space atlas labels
    ├── parc/
    │   └── <ent>_seg-<seg>_level-<N>_from-ABAv3_dseg.nii.gz   # subject-space labels (download)
    ├── featuremap/
    │   └── <ent>_stain-<stain>_space-ABAv3_level-<N>_fieldfrac.nii.gz
    └── qc/
        ├── <ent>_stain-<stain>_level-<N>_space-ABAv3_regqc.html
        └── <ent>_qc.json                               # scan QC metrics (from qc_scan)
```

Notes:

- The wide `mergedsegstats.tsv` keeps its `<stain>+<metric>` columns (primary methods
  only) for existing users and scripts. The explorer and group analysis read only
  `roistats.parquet`.
- `featuremap/` holds the template-space field-fraction maps, including the CD31 vessel
  field fraction. They are for download and for group maps; the
  explorer colours regions from the tables, not from these maps.
- **View copies** are OME-Zarr 0.5 (Zarr v3) stores, so the same files open in
  Neuroglancer, napari and zarrnii, in the coordinate frames of
  [Coordinate systems](#coordinate-systems). Their encoding is for viewing, not analysis:
  - intensities as uint8 on a square-root scale between the 0.5th and 99.99th
    percentiles (the percentiles are recorded in the store's attributes);
  - chunks thin in z for xy viewing (256×256×16, x×y×z), in shards of about 512 MB;
  - a multiscale pyramid down to a few hundred voxels across;
  - lossless `zstd` chunk compression, which stock Neuroglancer reads;
  - masks and labels as unsigned integers at full resolution with a nearest-neighbour
    pyramid; region names come from the release's label table.

  The chunk codec is a parameter of the view encoding: a lossy codec (JPEG or JPEG XL,
  which would need a matching Zarr codec in our Neuroglancer build) or a Neuroglancer
  precomputed copy can be added later if measurements show the explorer needs fewer
  bytes per screen.

### Release root and storage format

```
--release_root PATH | mind://...    default: <output_dir>/release
--release_zarr ome.zarr | ozx       default: ome.zarr
```

- **Release root.** The release can be written to a local path or straight to MIND. MIND's
  service jobs on ingested data write to a `mind://` root (the run's own derivative
  dataset), so there is no separate upload step. Only the release is written there: the
  run's other outputs (full-resolution masks, template-space images, work files) stay
  in `<output_dir>` and the release step reads them from there.
- **Zarr format.** View copies are sharded OME-Zarr directory stores by default, with the
  chunking and sharding of [Layout](#layout). `ozx` (zipped OME-Zarr, unsharded) is
  allowed for local releases only, e.g. on shared filesystems with inode limits; it is an
  error with a `mind://` root. A release bound for MIND is always sharded `ome.zarr`,
  because Neuroglancer reads directory stores by HTTP range requests. A local `ozx`
  release is converted to `ome.zarr` when it is uploaded (or refused by the upload).
  The format is recorded in the provenance sidecar (`ZarrFormat`).
- **Completion.** The provenance sidecar is written last; its `Files` list (sizes and
  sha256, computed while writing) is what marks the release complete. MIND registers the
  derivative dataset only when the sidecar exists and the files match it, so an
  interrupted job leaves an unregistered prefix to clean up, never a partial release.
- **Writing to `mind://`:**
  - each view-copy shard is written by exactly one task (dask chunks aligned to shards),
    and shards stay under the MIND proxy's maximum single upload (it has no multipart
    upload);
  - rules that write to `mind://` declare a small local receipt file (URI, sha256) as
    their Snakemake output, so completion is tracked without a storage plugin;
  - jobs refresh their MIND token while they run, since encoding a large subject can
    outlast one token.

**Where each release file comes from** (SPIMquant output names, paths relative to
`<output_dir>`; `<p>` = the stain's primary seg method, `<o>` = its primary object
method):

| Release file | SPIMquant output |
| --- | --- |
| `tabular/<ent>_space-ABAv3_roistats.parquet` | Melted from `tabular/<ent>_seg-<seg>_from-ABAv3_mergedsegstats.tsv` (primary methods); comparison methods from `tabular/<ent>_stain-<stain>_seg-<seg>_from-ABAv3_level-<N>_obj-<obj>_desc-<seg_method>_segstats.tsv` |
| `tabular/<ent>_seg-<seg>_from-ABAv3_mergedsegstats.tsv` | Same name |
| `tabular/<ent>_space-ABAv3_instances.parquet` | `tabular/<ent>_space-ABAv3_regionprops.parquet`, plus `label`, `instance_id` and the subject columns; `stain` renamed `channel` |
| `tabular/<ent>_space-ABAv3_coloc.parquet` | Same name |
| `micr/<ent>_SPIM.json` | The input BIDS sidecar of the scan |
| `view/<ent>_stain-<stain>_SPIM.ome.zarr` | Encoded from the input image |
| `view/<ent>_stain-<stain>_space-ABAv3_SPIM.ome.zarr` | `micr/<ent>_stain-<stain>_space-ABAv3_res-<r>um_level-<N>_desc-<desc>_SPIM.ozx` (`--template_spim_res`), re-encoded |
| `xfm/<ent>_from-<acq>_to-<ref acq>_xfm.txt` | `xfm/<ent>_from-ref_to-subject_desc-<method>_type-ras_xfm.txt` (`--aux_spim`) |
| `featuremap/<ent>_stain-<stain>_space-ABAv3_level-<N>_fieldfrac.nii.gz` | `seg/<ent>_stain-<stain>_space-ABAv3_level-<N>_desc-<p>_fieldfrac.nii.gz` (CD31: `vessels/..._desc-vesselfm_fieldfrac.nii.gz`) |
| `view/<ent>_stain-<stain>_mask.ome.zarr` | `seg/<ent>_stain-<stain>_level-<segmentation_level>_desc-<p>_mask.ozx`, re-encoded |
| `view/<ent>_seg-all_from-ABAv3_dseg.ome.zarr` | The template's `seg-all` labels warped to the scan at `segmentation_level` (today SPIMquant warps them to `registration_level` only) |
| `parc/<ent>_seg-<seg>_level-<N>_from-ABAv3_dseg.nii.gz` | Same name |

### Regional table (`roistats.parquet`)

The long form of `mergedsegstats.tsv`: one row per participant × parcellation × label ×
channel × method × metric. It holds the primary methods and any comparison methods
(as extra rows with their own `seg_method` / `object_method`).

| Column | Type | Description |
| --- | --- | --- |
| `participant_id` | string | `sub-<label>` |
| `session_id` | string, nullable | `ses-<label>` |
| `acq` | string, nullable | The scan's `acq-` entity |
| `scan_role` | string | `reference` or `additional` (see [Scans](#scans)) |
| `modality` | string | `micr` for SPIM rows (MRI rows use `anat`, `dwi`, `func`) |
| `channel` | string | Stain name, or `coloc` for colocalization metrics |
| `seg_method` | string | Mask method (e.g. `lantern`); null for coloc rows |
| `object_method` | string, nullable | `components` or `watershed`; null for mask-only metrics (`fieldfrac`, `volume`) |
| `primary` | bool | True for the stain's primary methods |
| `template` | string | `ABAv3` |
| `seg` | string | Parcellation (`all`, `coarse`, `mid`, `fine`) |
| `label` | int32 | `index` in `tpl-ABAv3_seg-<seg>_dseg.tsv` |
| `metric` | string | Metric name from the vocabulary below |
| `value` | float64 | Metric value; null if not measurable (e.g. region outside field of view) |
| `units` | string | Units of `value` |

Region names are not repeated in every row. Readers join `label` to the release's
label table.

**Metric vocabulary (SPIM, schema 0.1).** These are the metrics SPIMquant writes today
(`seg_metrics`, `coloc_seg_metrics`, and the per-instance features averaged by
`merge_into_segstats_tsv.py`):

| `metric` | Units | Meaning | Source |
| --- | --- | --- | --- |
| `fieldfrac` | percent | Mean percentage of the region occupied by the mask | `fieldfracstats` |
| `count` | objects | Objects whose centroid falls in the region | `countstats` |
| `volume` | mm³ | Region volume in subject space | `fieldfracstats` |
| `density` | objects/mm³ | `count / volume` | computed |
| `nvoxels` | voxels | Total voxels of objects in the region | `regionpropstats` (sum) |
| `<stain>_sdt` | mm | Mean signed distance from object centroids to the vessel wall (negative = inside), vessel stain only | instance table (mean) |
| `<stain>_dist` | mm | Mean centroid distance to the nearest object of `<stain>` | instance table (mean) |
| `<stain>_edgedist` | mm | Mean approximate distance to the surface of the nearest object of `<stain>` | instance table (mean) |
| `objvolume` | mm³ | Mean object volume (the objects' `volume`, renamed per region) | instance table (mean) |
| `mean_intensity`, `max_intensity`, `peak_intensity`, `elongation` | a.u. / ratio | Mean over the region's objects (`watershed` only) | instance table (mean) |
| `coverage` | fraction | Fraction of the region (as labelled in the reference scan) inside the field of view; additional scans only | `coveragestats` |
| `overlapratio`, `distance` | ratio / mm | Colocalization metrics (`channel = coloc`) | `colocsegstats` |

A metric not in this list is allowed only with a description in the sidecar and a minor
schema bump.

### Instance table (`instances.parquet`)

One row per detected object (plaque, microglial cell, ...) for every stain with
objects, using its primary methods. SPIMquant writes most of it as
`tabular/<ent>_space-ABAv3_regionprops.parquet`: the primary objects of every stain,
with `stain`, `seg_method`, `object_method`, `volume` and the distance columns. The
release adds the ABAv3 `seg-all` `label`, `instance_id` and the subject columns, and
renames `stain` to `channel`.

| Column | Type | Description |
| --- | --- | --- |
| `participant_id` | string | `sub-<label>` |
| `session_id` | string, nullable | `ses-<label>` |
| `acq` | string, nullable | The scan's `acq-` entity |
| `scan_role` | string | `reference` or `additional` (see [Scans](#scans)) |
| `channel` | string | Stain (SPIMquant's `stain` column) |
| `seg_method` | string | Mask method |
| `object_method` | string | `components` or `watershed` |
| `instance_id` | int64 | Unique within (`participant_id`, `session_id`, `acq`, `channel`); stable within a release |
| `label` | int32 | ABAv3 `seg-all` index of the region containing the centroid (0 = outside the brain) |
| `template_x`, `template_y`, `template_z` | float64, mm | Centroid in ABAv3 space |
| `pos_x`, `pos_y`, `pos_z` | float64, mm | Centroid in subject physical space |
| `nvoxels` | float64 | Object size in voxels at `segmentation_level` |
| `volume` | float64, mm³ | `nvoxels` × voxel volume, so sizes compare across scans |
| `<stain>_dist` | float64, mm | One per stain with objects: centroid distance to the nearest object of that stain (nearest other object for its own stain) |
| `<stain>_edgedist` | float64, mm | One per stain with objects: `<stain>_dist` minus the nearest object's equivalent-sphere radius, clipped at 0 |
| `<stain>_sdt` | float64, mm | Vessel stain only: signed distance to the vessel wall, negative inside |
| `seed_x`, `seed_y`, `seed_z` | float64, mm, nullable | Seed position (`watershed` only) |
| `mean_intensity`, `max_intensity`, `peak_intensity` | float64, nullable | Intensities of the bias-corrected image (`watershed` only) |
| `elongation` | float64, nullable | 1 for a sphere (`watershed` only) |

- The `label` is looked up in subject space, from the subject-space `seg-all` labels and
  `pos_*`, the same way `map_atlas_to_regionprops.py` labels objects for the regional
  tables. Instance counts per label therefore match `count` in `roistats.parquet`.
- Rows are sorted by `channel`, `label`. Row groups hold at most 1,000,000
  rows, so Parquet statistics let DuckDB-WASM skip row groups when filtering by stain or
  region.
- Distance column names use the stain names of the subject, so they differ between
  datasets; the sidecar lists them (`"DistanceStains"`, `"VesselStain"`).

### Colocalization table

`coloc.parquet` is released as SPIMquant writes it (`object_id_a`, `object_id_b`,
`stain_a`, `stain_b`, `distance`, `overlapratio`, `template_coloc_*`, ...), plus
`participant_id` and `session_id`. It pairs every two stains with objects, whatever
their methods, and exists when a subject has at least two.

### Vessel graph (planned, not in 0.1)

The CD31 skeleton graph (`nodes.parquet`, `edges.parquet`) becomes the CD31 instance
table in a later minor version: one row per vessel segment (edge), with template
coordinates of its endpoints, length, calibre (radius), and `seg-all` label. Until then
the graph stays out of the release.

### Provenance sidecar (`<ent>_provenance.json`)

One per participant release (`sub-<label>[_ses-<label>]_provenance.json`), covering all
its scans. It is what PR2's compatibility check reads, and what MIND records for the
derivative.

```json
{
  "SchemaVersion": "0.1.0",
  "ReleaseType": "participant",
  "Created": "2026-11-21T14:00:00-05:00",
  "Pipeline": {
    "Name": "spimquant",
    "Version": "0.11.0",
    "CodeURL": "https://github.com/khanlab/SPIMquant/tree/<commit>",
    "Commit": "<sha>",
    "Dirty": false,
    "SnakebidsRunId": "20261106T134023-25f2"
  },
  "ConfigHash": "sha256:<hex>",
  "ZarrFormat": "ome.zarr",
  "Template": {
    "Name": "ABAv3",
    "Source": "https://zenodo.org/records/18906782",
    "LabelTables": {"all": "sha256:<hex>", "coarse": "sha256:<hex>"}
  },
  "Registration": {"Method": "<desc of the selected warp>", "Level": 5, "Stain": "PI",
                   "Orientation": "<--orientation>", "TemplateCrop": null},
  "Reused": [
    {"Outputs": ["brainmask", "transforms"],
     "From": "mind://<dataset>/derivatives/spimquant_reg",
     "Pipeline": {"Name": "spimquant", "Version": "0.10.2", "Commit": "<sha>",
                  "SnakebidsRunId": "<run id of the prior dataset>"},
     "ConfigHash": "sha256:<hex>",
     "RegQC": "qc/<ent>_stain-<stain>_level-<N>_space-ABAv3_regqc.html"}
  ],
  "Segmentation": [
    {"Stain": "Abeta", "SegMethod": "lantern", "ObjectMethod": "components", "Primary": true,
     "Source": "default", "Model": "<lantern_config>", "Level": 0},
    {"Stain": "Iba1", "SegMethod": "gmm+n3k1", "ObjectMethod": "watershed", "Primary": true,
     "Source": "default", "Level": 0}
  ],
  "Vessels": {"Stain": "CD31", "Method": "vesselfm", "Level": 0},
  "Scans": [
    {"Role": "reference",
     "Entities": {"subject": "01", "session": "01", "sample": "brain", "acquisition": "1x"},
     "Stains": ["Abeta", "Iba1", "CD31", "PI"],
     "URI": "mind://<dataset>/sub-01/ses-01/micr/...acq-1x_SPIM.ome.zarr",
     "ContentFingerprint": null},
    {"Role": "additional", "Reference": "1x",
     "Entities": {"subject": "01", "session": "01", "sample": "brain", "acquisition": "12x"},
     "Stains": ["Abeta", "Iba1", "CD31", "PI"],
     "URI": "mind://<dataset>/sub-01/ses-01/micr/...acq-12x_SPIM.ome.zarr",
     "ContentFingerprint": null,
     "RegistrationToReference": {"Method": "affine", "Init": "match", "Stains": ["PI"],
                                 "QC": "qc/<ent>_from-ref_to-subject_regqc.html"}}
  ],
  "Files": [
    {"path": "tabular/sub-01_ses-01_space-ABAv3_roistats.parquet", "bytes": 123456, "sha256": "<hex>"}
  ]
}
```

| Field | Rule |
| --- | --- |
| `SchemaVersion` | The schema version of this spec the release was written for. |
| `Pipeline` | Taken from snakebids' `code/snakebids/runs/<run_id>/provenance.json` of the run that produced the outputs. |
| `ConfigHash` | sha256 of the resolved config as canonical JSON (sorted keys), after removing keys that do not change results: input and output paths, `work`, dask/cluster and resource settings, `--until`/targets. The excluded keys are listed in the sidecar's JSON schema, not per file. |
| `Template` | Name, source record, and sha256 of each released label table. |
| `Segmentation` | One entry per (stain, seg method, object method) in the tables, with whether it is primary and whether it came from the CLI or a default. |
| `Scans` | One entry per scan in the release: its role, entities, stains and input. An additional scan names its reference scan's `acq` and its registration to it (`--aux_reg_method`, `--aux_reg_init`, `--aux_reg_stains`) with the QC report. |
| `Scans[].URI` | The input path or `mind://` URI as given to SPIMquant. |
| `Scans[].ContentFingerprint` | Left null by SPIMquant. The fingerprint of the level-0 data is computed by `mind_client` (from the OME-Zarr or the `.ims`) and filled in on upload. |
| `Registration` | How the subject was registered to the template, including the settings that change the result without changing file names (`--orientation`, `--template_crop`, `--sloppy`). |
| `Reused` | Empty unless outputs came from a prior dataset (`--reuse_derivatives`). One entry per prior dataset: which outputs, where from, and the prior run's pipeline version, commit, run id and config hash, read from that dataset's `code/snakebids/runs/`. `Pipeline` and `ConfigHash` above describe the run that made the release. |
| `Files` | Every file in the participant release except this sidecar, with size and sha256, so the upload can be verified. The sidecar is written last and marks the release complete (see [Release root and storage format](#release-root-and-storage-format)). |
| `ZarrFormat` | `ome.zarr` or `ozx`: how the view copies are stored (`--release_zarr`). |

### Reused outputs (`--reuse_derivatives`)

Reusing outputs (e.g. brain masks and registrations) from a prior SPIMquant dataset is a
way to split work between runs; it does not change what a release is. Reused outputs
are not shared or referenced between releases: each participant release holds its own
copy, so it stays self-contained. The run that makes the release:

- writes no symlinks into a release: reused files (also those linked with
  `--reuse_symlink` in the run's output) are copied in as real files, so the release is
  the complete upload set;
- includes the prior dataset's registration QC report, since SPIMquant does not make
  one for reused transforms, so every participant release has one;
- records the prior run in `Reused`, so the provenance names every run whose outputs
  are in the release;
- checks the settings SPIMquant does not check when reusing (`--orientation`,
  `--template_crop`, `--sloppy`) against the prior run's config, and refuses to write a
  release if they differ.

Release consumers treat a release with reused outputs like any other; the compatibility
check of the aggregation manifest compares the registration settings in `Registration`,
which are the same whether the registration was reused or not.

---

## Group release

Written by the group level under `release/group/`. It aggregates participant releases
(local or `mind://`) listed in the aggregation manifest, or found under `root` when no
manifest is given. Its inputs are usually the ingestion-time participant releases in
MIND, possibly from several datasets; it adds the covariates and QC status that
participant releases leave out.

```
release/
├── participants.tsv                       # participant_id, session_id + covariates
├── participants.json                      # covariate descriptions (BIDS style)
├── release_manifest.json                  # explorer manifest (below)
└── group/
    ├── space-ABAv3_roistats.parquet       # union of participant roistats
    ├── space-ABAv3_instances.parquet      # union of participant instance tables
    ├── space-ABAv3_coloc.parquet
    ├── seg-<seg>_from-ABAv3_allsubjects.tsv               # wide, for users
    ├── seg-<seg>_from-ABAv3_contrast-<c>_groupstats.tsv   # when contrasts are requested
    └── maps/
        ├── seg-<seg>_space-ABAv3_metric-<stain>+<metric>_mean.nii.gz
        └── seg-<seg>_space-ABAv3_contrast-<c>_metric-<stain>+<metric>_<stat>.nii.gz
```

- Group tables are plain unions of participant tables, all scans included; they carry
  no covariates. Consumers filter on `acq` / `scan_role` and never pool a subject's
  scans; group statistics and maps are made per acquisition, as SPIMquant does today
  (reference scans unprefixed, additional scans named with their `acq-`).
  Covariates live once in `participants.tsv` and are joined by `participant_id` and
  `session_id`.
- `instances.parquet` is sorted by `channel`, `label`, `participant_id`, with the
  same row-group limit as the participant table.
- Group maps are the template-space region maps SPIMquant already paints
  (`featuremap`, group-stats maps), plus a per-region group mean, for primary methods.
- Participant releases are not copied into the group release. The manifest points to
  them.

### Release manifest (`release_manifest.json`)

The entry point for the explorer and for DOIs. It lists what the release holds; the
explorer reads it first and then reads tables by range.

```json
{
  "SchemaVersion": "0.1.0",
  "ReleaseId": "<dataset>-2026-11-21",
  "Title": "...",
  "Created": "2026-11-21T14:00:00-05:00",
  "DOI": null,
  "License": "CC-BY-4.0",
  "Template": {"Name": "ABAv3", "Source": "https://zenodo.org/records/18906782"},
  "LabelTables": {"all": "tpl-ABAv3/tpl-ABAv3_seg-all_dseg.tsv", "coarse": "..."},
  "Subjects": [
    {"participant_id": "sub-01", "session_id": "ses-01", "MRI": ["T2w", "dwi"],
     "Release": "sub-01/ses-01/",
     "Scans": [
       {"acq": "1x", "Role": "reference", "Stains": ["Abeta", "Iba1", "CD31"], "QC": "pass",
        "Microscope": "UltraMicroscope Blaze II", "Magnification": 1.1},
       {"acq": "12x", "Role": "additional", "Reference": "1x", "Coverage": 0.08,
        "Stains": ["Abeta", "Iba1", "CD31"], "QC": "pass",
        "Microscope": "UltraMicroscope Blaze II", "Magnification": 12}
     ]}
  ],
  "Stains": [
    {"Name": "Abeta", "SegMethod": "lantern", "ObjectMethod": "components", "Target": "amyloid plaques",
     "PrimaryAntibody": "<RRID or vendor, catalog no.>", "SecondaryAntibody": "<...>",
     "Fluorophore": "<...>"},
    {"Name": "Iba1", "SegMethod": "gmm+n3k1", "ObjectMethod": "watershed", "Target": "microglia"},
    {"Name": "CD31", "SegMethod": "vesselfm", "ObjectMethod": null, "Target": "vessels"}
  ],
  "Metrics": [{"Name": "density", "Units": "objects/mm^3", "Description": "..."}],
  "MRI": {
    "Contrasts": ["T2w", "dwi"],
    "Regional": [{"Pipeline": "<name>", "Version": "<v>", "Metrics": ["fa", "md"]}],
    "Connectivity": [{"Pipeline": "<name>", "Version": "<v>", "Type": "structural",
                      "Measure": "streamlinecount", "Seg": "coarse"}]
  },
  "Pipelines": [{"Name": "spimquant", "Version": "0.11.0", "Subjects": 24}],
  "Tables": {
    "roistats": "group/space-ABAv3_roistats.parquet",
    "instances": "group/space-ABAv3_instances.parquet",
    "mri_roistats": "group/space-ABAv3_mri_roistats.parquet",
    "connectivity": ["group/seg-coarse_space-ABAv3_meas-streamlinecount_relmat.parquet"]
  },
  "Files": [{"path": "group/space-ABAv3_roistats.parquet", "bytes": 0, "sha256": "<hex>"}]
}
```

- `DOI` is null until one is minted; minting does not change any other file.
- `Subjects[].Scans[].QC` comes from the MIND catalog (QC status per acquisition and
  run). An explorer release includes `pass` scans only.
- `Scans[].Coverage` is the fraction of the brain inside an additional scan's field of
  view.
- `Stains[]` antibody, fluorophore and `Target` values, and `Scans[]` microscope and
  magnification, are summarised from the scans' `SPIM.json` sidecars (see
  [Sample and acquisition metadata](#sample-and-acquisition-metadata)). A stain whose
  scans disagree lists every value.

---

## MRI derivative schema

Any pipeline (ours or a user's) can produce MRI derivatives that MIND accepts, as long as
they follow this layout. They are a BIDS-derivative dataset per pipeline:

```
<pipeline>/
├── dataset_description.json               # GeneratedBy: name, version, CodeURL
└── sub-<label>/[ses-<label>/]<datatype>/  # anat, dwi, func
    ├── <ent>_space-ABAv3_roistats.parquet
    ├── <ent>_space-ABAv3_roistats.json
    ├── <ent>_seg-<seg>_space-ABAv3_meas-<measure>_relmat.parquet
    ├── <ent>_seg-<seg>_space-ABAv3_meas-<measure>_relmat.json
    └── <ent>_provenance.json
```

### Regional table

The same columns as the SPIM [regional table](#regional-table-roistatsparquet), with:

- `modality` = the BIDS datatype (`anat`, `dwi`, `func`);
- `channel` = the contrast or BIDS suffix the metric was measured on (`T2w`, `dwi`,
  `bold`);
- `seg_method` = the pipeline's own label for the method (may be null), `object_method`
  null, `primary` true;
- `label` = index in the ABAv3 label table for `seg`.

Every metric used must be described in the sidecar with `LongName`, `Description` and
`Units`. A metric vocabulary for MRI is not fixed in 0.1; pipelines add entries as they
mature.

### Connectivity (`relmat.parquet`)

An edge list. Naming follows the BIDS connectivity proposal (BEP017): suffix `relmat`,
entity `meas-` for the edge measure.

| Column | Type | Description |
| --- | --- | --- |
| `participant_id` | string | `sub-<label>` |
| `session_id` | string, nullable | `ses-<label>` |
| `src_label` | int32 | Node index in the `seg` label table |
| `dst_label` | int32 | Node index in the `seg` label table |
| `weight` | float64 | Edge weight |

Sidecar (required fields):

| Field | Values |
| --- | --- |
| `ConnectivityType` | `structural`, `functional`, `effective` |
| `EdgeWeight` | What `weight` is, e.g. `streamline count`, `mean FA along streamlines`, `Pearson r`, `Fisher z` |
| `Units` | Units of `weight`, or `none` |
| `Directed` | `false`: each pair stored once with `src_label < dst_label`; `true`: both directions stored |
| `SelfConnections` | Whether rows with `src_label = dst_label` are present |
| `MissingEdge` | `zero` (absent rows mean weight 0) or `null` (absent rows mean not computed) |
| `Template`, `Seg` | `ABAv3` and the parcellation; nodes are the labels of `tpl-ABAv3_seg-<seg>_dseg.tsv` |
| `Nodes` | Optional list of label indices used, when the pipeline uses a subset of the table |

Nodes are ABAv3 labels, either the full `seg-all` set or a coarser declared parcellation
(`coarse`, `mid`, `fine`, `roi22`, `roi82`, `roi198`). Using the same nodes as SPIM
lets the explorer plot pathology at a node against that node's connectivity.

### MRI provenance

`<ent>_provenance.json` with the same fields as the SPIM sidecar where they apply:
`SchemaVersion`, `Pipeline` (name, version, code URL), `ConfigHash` (optional),
`Template`, the parcellations used (`Seg` with label-table sha256), `Input` (entities,
URI), and `Files`.

---

## Aggregation manifest (`--group-manifest`)

A TSV with one row per subject and session that group analysis should include. It
replaces finding subjects under `root` and reading covariates from
`bids_dir/participants.tsv`. Without `--group-manifest`, group-level behaviour does not
change.

```tsv
participant_id	session_id	derivatives	genotype	age_weeks	sex
sub-01	ses-01	/scratch/runA/release/sub-01/ses-01	APP	24	F
sub-02	ses-01	mind://projectX/derivatives/spimquant/sub-02/ses-01	WT	24	M
sub-03		../runB/release/sub-03	APP	30	F
```

| Column | Required | Rule |
| --- | --- | --- |
| `participant_id` | yes | `sub-<label>` |
| `session_id` | yes (may be empty) | `ses-<label>`, or empty for datasets without sessions |
| `derivatives` | yes | Path to a participant release directory (absolute, or relative to the manifest file), or a `mind://` URI to one |
| any other column | no | Covariate. Copied to the group `participants.tsv`; usable in `--group-stats-model`, `--group-stats-pairwise`, `--group-stats-within`, `--group-stats-where` |

- An optional `<manifest>.json` describes covariate columns (BIDS `participants.json`
  style: `Description`, `Levels`, `Units`).
- `(participant_id, session_id)` must be unique. Every scan in the referenced
  participant release is used; group outputs are made per acquisition (see
  [Scans](#scans)).
- `mind://` rows are staged by fetching the release tables and provenance sidecar only
  (a few MB per subject), not images.
- The first group rule checks compatibility from the provenance sidecars of all rows:
  - **fail** if `Template.Name`, any label-table sha256, the primary
    (stain, seg method, object method) entries, or the major `SchemaVersion` differ
    between rows;
  - **warn** if `Pipeline.Version` or `ConfigHash` differ, or if `Registration`
    differs (method, level or orientation), including through `Reused`.
- A participant release missing a table that the requested group outputs need is an
  error naming the row.

---

## Validation

Besides the schema checks below, the CI check samples objects from each instance table
and checks that their coordinates land on the expected label and mask values in every
view copy of the same space (template and subject), so a coordinate-frame mismatch
between OME-Zarr, NIfTI and the tables fails the build.

PR1 adds a JSON Schema for the provenance sidecar and the release manifest, and a
column-schema check for each table (names, types, units, allowed `metric` values). CI
runs the `release` target on `tests/bids_ds` and validates the output. MIND runs the same
checks on upload, for SPIM and MRI derivatives alike.

---

## Changes needed in SPIMquant

**Done in #248** (per-stain methods):

- Per-stain `--seg_method` / `--object_method` with primaries; `desc-`/`obj-` only on
  per-stain files; `cc`, `inst` and `+inst` rejected; GPU methods fall back without a
  GPU (`--gpu`, `gpu_seg_method_fallbacks`).
- Colocalization and the merged tables combine stains whatever their methods; fibre
  ROI metrics are in `mergedsegstats.tsv`.
- `<stain>_dist`, `<stain>_edgedist` for every stain with objects; `<vessel>_sdt`
  replaces `vesseldist` (mm).
- Objects carry `stain`, `seg_method`, `object_method` and `volume` (mm³); per region,
  `objvolume`.
- Column descriptions give `pos_*`, distances and colocalization radii in mm;
  `outputs.md` describes `regionpropstats.tsv` as one row per object.

**Still needed for the `release` target (PR1):**

| Item | Today (after #248) | Needed |
| --- | --- | --- |
| Long regional table | Wide `mergedsegstats.tsv` per `seg`; comparison methods in per-stain `segstats.tsv` | Write `roistats.parquet` (all segs and methods, with `primary`) |
| Instance label | `space-ABAv3_regionprops.parquet` has no atlas label; per-object labels only in the per-stain `regionpropstats.tsv` | Add `seg-all` `label`, using the same subject-space lookup |
| Instance identity | No id (coloc builds `object_id` from row order) | Add `instance_id`, `participant_id`, `session_id`; rename `stain` to `channel` |
| Provenance sidecar | Run provenance in `code/snakebids/runs/`; method plan printed at start-up | Write `<ent>_provenance.json` per participant, including the method plan |
| View copies | Template-space OME-Zarr at `--template_spim_res`, not encoded for viewing; no subject-space copies | OME-Zarr view copies (subject and template space, masks, labels) with the view encoding and the frames of [Coordinate systems](#coordinate-systems) |
| Scans | `--aux_spim` outputs named with `acq-`; group outputs per acquisition | `acq` and `scan_role` columns; `Scans` in the provenance sidecar; the additional scan's transform as a 4×4 affine JSON |
| Reused outputs | `--reuse_symlink` links files; no registration QC report for reused transforms; prior run not recorded | Real files in the release, the prior dataset's QC report copied, `Reused` in provenance, unchecked settings compared |
| Sample and acquisition metadata | In the input BIDS sidecars, not carried to outputs | Copy each scan's `SPIM.json`; validate the per-channel lists. SPIMpack adds the extension fields (`SampleFluorophore`, `ExcitationWavelength`, `EmissionWavelength`, `LaserPower`, `ExposureTime`) |
| Release tree and manifest | None | `release/` layout, label tables, `release_manifest.json` |
| Release root and format | Outputs under `<output_dir>` only; zarr outputs `ozx` by default (`zarr_store_ext`), sharding off for `ozx` | `--release_root` (local or `mind://`) and `--release_zarr` (`ome.zarr` default, `ozx` local only); receipt files for `mind://` outputs; provenance sidecar written last |

---

## Open questions

1. **Watershed cost.** `components` is the default; whether `watershed` should be
   the default for cell stains such as Iba1 waits on a runtime comparison on one
   subject at level 0.
2. **Group maps.** Are per-region group means the only maps needed, or also voxel-wise
   group averages of the template-space field-fraction and count maps?
3. **Template-space view copy resolution.** Which `--template_spim_res` the explorer
   uses (one resolution, or a pyramid).
4. **MRI regional metrics and connectivity types that go in first** (also open in the
   portal plan).
5. **Lossy view encoding.** Whether lossless uint8 `zstd` is small enough for the
   explorer over the network, or a JPEG / JPEG XL Zarr codec in our Neuroglancer build
   (or a precomputed copy) is needed. Measure bytes and load time per screen against the
   existing precomputed copy of `test_mind_4x`.
6. **Shard size and the MIND proxy.** The ~512 MB view-copy shards must fit the proxy's
   (and the backing object store's) maximum single upload; check the limit before
   fixing the shard size.
7. **OME-Zarr transforms.** When the OME-Zarr coordinate-transformations proposal lands,
   the additional scan's affine can move from a sidecar JSON into the store.
