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

## Conventions

| Item | Rule |
| --- | --- |
| Template | `ABAv3` for every release that MIND accepts. Template-space data and template coordinates are in ABAv3 physical space. |
| Label table | `tpl-ABAv3_seg-<seg>_dseg.tsv` from the template's Zenodo record (columns `index`, `name`, ...). The release copies the label tables it uses, so a release is self-contained. |
| Canonical labels | `seg-all` (the full ABAv3 label set). Instance tables carry `seg-all` labels only. Coarser parcellations come from the label hierarchy or from the regional table. |
| Released parcellations | Regional tables cover the template's `default_segs` (`all`, `coarse`, `mid`, `fine` for ABAv3). |
| Units | Coordinates and distances in **mm**; volumes in **mm³**; densities in **objects/mm³**; field fraction in **percent (0–100)**; intensities in arbitrary units. Units are written in the table (`units` column) or in the JSON sidecar, never assumed. |
| Subject key | `participant_id` (`sub-<label>`) and `session_id` (`ses-<label>`, or null when the dataset has no sessions). Every table carries both, so tables from many subjects concatenate without joins. |
| Other BIDS entities | `sample`, `acq` and similar appear in filenames and in the provenance sidecar. One release holds one acquisition per subject and session. |
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
| Iba1 | Microglia: `gmm+n3k1` masks by default, split by `watershed` so touching cells count separately | Regional table rows; instance table rows (one per cell); `Iba1_dist` and `Iba1_edgedist` on every object |
| CD31 | Vessels: `vesselfm` mask and its signed distance transform (SDT) | `CD31_sdt` on every Abeta and Iba1 object, and its regional mean; template-space vessel field-fraction map. The vessel graph (with calibre per edge) becomes the CD31 instance table later. |

Abeta and Iba1 objects are colocalized with each other even though they use different
methods.

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
    built-in defaults: Abeta=lantern (also abeta, BetaAmyloid); every other stain gmm+n3k1

--object_method SPEC [SPEC ...]
    STAIN=OBJ[,OBJ...] | OBJ   same rules as --seg_method
    OBJ: watershed | components | none
    built-in defaults: Abeta=components; every other stain watershed

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
- **File names:** per-stain files carry the method (`desc-<seg_method>` on masks and
  mask-derived maps; `desc-<seg_method>_obj-<object_method>` on per-stain object
  tables). Cross-stain results and release files carry neither; the methods are in
  table columns and sidecars.
- Vessels keep their own options (`--stains_for_vessels`, `--no_vessels`) and outputs
  (mask, SDT, graph); they have no object method.

---

## Participant release

### Layout

The release is a BIDS-derivative dataset of its own, written under `<output_dir>/release/`.
Paths below are relative to that root; `<ent>` is the subject's entity prefix
(`sub-01_ses-01_sample-brain`, without `stain-`).

```
release/
├── dataset_description.json              # DatasetType: derivative, GeneratedBy: spimquant
├── tpl-ABAv3/
│   ├── tpl-ABAv3_seg-<seg>_dseg.tsv      # label tables used by this release
│   └── tpl-ABAv3_seg-<seg>_dseg.json
└── sub-<label>/[ses-<label>/]
    ├── <ent>_provenance.json             # provenance sidecar (below)
    ├── tabular/
    │   ├── <ent>_space-ABAv3_roistats.parquet          # regional table, long (all segs, stains, descs)
    │   ├── <ent>_space-ABAv3_roistats.json
    │   ├── <ent>_seg-<seg>_from-ABAv3_mergedsegstats.tsv   # wide TSV, primary methods, for users
    │   ├── <ent>_seg-<seg>_from-ABAv3_mergedsegstats.json
    │   ├── <ent>_space-ABAv3_instances.parquet         # instance table (all stains, primary methods)
    │   ├── <ent>_space-ABAv3_instances.json
    │   ├── <ent>_space-ABAv3_coloc.parquet             # colocalization pairs (2+ stains with objects)
    │   └── <ent>_space-ABAv3_coloc.json
    ├── view/                                           # Neuroglancer view copies
    │   ├── <ent>_stain-<stain>_SPIM.precomputed/        # subject space, full resolution pyramid
    │   ├── <ent>_stain-<stain>_space-ABAv3_SPIM.precomputed/
    │   ├── <ent>_stain-<stain>_mask.precomputed/        # subject-space primary masks
    │   └── <ent>_seg-all_from-ABAv3_dseg.precomputed/   # subject-space atlas labels
    ├── parc/
    │   └── <ent>_seg-<seg>_level-<N>_from-ABAv3_dseg.nii.gz   # subject-space labels (download)
    ├── featuremap/
    │   └── <ent>_stain-<stain>_level-<N>_space-ABAv3_fieldfrac.nii.gz
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
- View copies are Neuroglancer precomputed, sharded. Intensity images are JPEG-encoded
  uint8 with square-root intensity scaling; masks and labels use
  `compressed_segmentation`. The label view carries segment properties built from the
  label table, so names show in the viewer. Encoding details follow `mind-ng viewcopy`
  (mind-neuroglancer, see its FINDINGS.md) and are fixed in PR1.

### Regional table (`roistats.parquet`)

The long form of `mergedsegstats.tsv`: one row per participant × parcellation × label ×
channel × method × metric. It holds the primary methods and any comparison methods
(as extra rows with their own `seg_method` / `object_method`).

| Column | Type | Description |
| --- | --- | --- |
| `participant_id` | string | `sub-<label>` |
| `session_id` | string, nullable | `ses-<label>` |
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
| `mean_intensity`, `max_intensity`, `peak_intensity`, `elongation` | a.u. / ratio | Mean over the region's objects (`watershed` only) | instance table (mean) |
| `coverage` | fraction | Fraction of the region inside the field of view (additional scans only) | `coveragestats` |
| `overlapratio`, `distance` | ratio / mm | Colocalization metrics (`channel = coloc`) | `colocsegstats` |

A metric not in this list is allowed only with a description in the sidecar and a minor
schema bump.

### Instance table (`instances.parquet`)

One row per detected object (plaque, microglial cell, ...) for every stain with
objects, using its primary methods. It is built from the aggregated template-space
`regionprops.parquet`, with distances and an ABAv3 `seg-all` label added downstream.

| Column | Type | Description |
| --- | --- | --- |
| `participant_id` | string | `sub-<label>` |
| `session_id` | string, nullable | `ses-<label>` |
| `channel` | string | Stain (today's `stain` column) |
| `seg_method` | string | Mask method |
| `object_method` | string | `components` or `watershed` |
| `instance_id` | int64 | Unique within (`participant_id`, `session_id`, `channel`); stable within a release |
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

One per participant release. It is what PR2's compatibility check reads, and what MIND
records for the derivative.

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
  "Template": {
    "Name": "ABAv3",
    "Source": "https://zenodo.org/records/18906782",
    "LabelTables": {"all": "sha256:<hex>", "coarse": "sha256:<hex>"}
  },
  "Registration": {"Method": "<desc of the selected warp>", "Level": 5, "Stain": "PI"},
  "Segmentation": [
    {"Stain": "Abeta", "SegMethod": "lantern", "ObjectMethod": "components", "Primary": true,
     "Source": "default", "Model": "<lantern_config>", "Level": 0},
    {"Stain": "Iba1", "SegMethod": "gmm+n3k1", "ObjectMethod": "watershed", "Primary": true,
     "Source": "default", "Level": 0}
  ],
  "Vessels": {"Stain": "CD31", "Method": "vesselfm", "Level": 0},
  "Input": {
    "Entities": {"subject": "01", "session": "01", "sample": "brain"},
    "Stains": ["Abeta", "Iba1", "CD31", "PI"],
    "URI": "mind://<dataset>/sub-01/ses-01/micr/...ome.zarr",
    "ContentFingerprint": null
  },
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
| `Input.URI` | The input path or `mind://` URI as given to SPIMquant. |
| `Input.ContentFingerprint` | Left null by SPIMquant. The fingerprint of the level-0 data is computed by `mind_client` (from the OME-Zarr or the `.ims`) and filled in on upload. |
| `Files` | Every file in the participant release except this sidecar, with size and sha256, so the upload can be verified. |

---

## Group release

Written by the group level under `release/group/`. It aggregates participant releases
(local or `mind://`) listed in the aggregation manifest, or found under `root` when no
manifest is given.

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

- Group tables are plain unions of participant tables; they carry no covariates.
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
    {"participant_id": "sub-01", "session_id": "ses-01",
     "Stains": ["Abeta", "Iba1", "CD31"], "MRI": ["T2w", "dwi"],
     "QC": "pass", "Release": "sub-01/ses-01/"}
  ],
  "Stains": [
    {"Name": "Abeta", "SegMethod": "lantern", "ObjectMethod": "components", "Target": "amyloid plaques"},
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
- `Subjects[].QC` comes from the MIND catalog (QC status per acquisition and run). An
  explorer release includes `pass` subjects only.
- `Stains[].Target` is free text supplied by whoever builds the release.

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
- `(participant_id, session_id)` must be unique.
- `mind://` rows are staged by fetching the release tables and provenance sidecar only
  (a few MB per subject), not images.
- The first group rule checks compatibility from the provenance sidecars of all rows:
  - **fail** if `Template.Name`, any label-table sha256, the primary
    (stain, seg method, object method) entries, or the major `SchemaVersion` differ
    between rows;
  - **warn** if `Pipeline.Version` or `ConfigHash` differ.
- A participant release missing a table that the requested group outputs need is an
  error naming the row.

---

## Validation

PR1 adds a JSON Schema for the provenance sidecar and the release manifest, and a
column-schema check for each table (names, types, units, allowed `metric` values). CI
runs the `release` target on `tests/bids_ds` and validates the output. MIND runs the same
checks on upload, for SPIM and MRI derivatives alike.

---

## Changes needed in SPIMquant

What PR1 has to add or fix for current outputs to meet this spec:

| Item | Today | Needed |
| --- | --- | --- |
| Methods per stain | `desc-<method>` keys every result; colocalization and merged tables only combine stains segmented by the same method; LANTERN special-cased to one stain | Per-stain `--seg_method` / `--object_method` with primaries; `desc-`/`obj-` only on per-stain files; `+inst` rejected (separate PR, before the `release` target) |
| Distances | `vesseldist` only (SDT of the vessel mask) | `<stain>_dist`, `<stain>_edgedist` for every stain with objects; `<vessel>_sdt` replaces `vesseldist` |
| Long regional table | Wide `mergedsegstats.tsv` per `seg` and method | Write `roistats.parquet` (all segs and methods in one file) |
| Instance label | Aggregated `regionprops.parquet` has no atlas label; per-object labels exist only in the per-stain `regionpropstats.tsv` | Add `seg-all` `label` to the instance table, using the same subject-space lookup |
| Instance identity | No id in `regionprops.parquet` (coloc builds `object_id` from row order) | Add `instance_id`, `participant_id`, `session_id`, `seg_method`, `object_method`, `volume` |
| `pos_*` units | Column descriptions say voxels; values are physical (they are warped with the physical-space transform) | Fix the descriptions to mm |
| `vesseldist` units | Described as voxels; the SDT is computed in physical units | Describe as mm (as `<vessel>_sdt`) |
| Provenance sidecar | Run provenance in `code/snakebids/runs/` | Write `<ent>_provenance.json` per participant |
| View copies | None | Neuroglancer precomputed copies (see layout) |
| Release tree and manifest | None | `release/` layout, label tables, `release_manifest.json` |
| Output docs | `outputs.md` calls `regionpropstats.tsv` one row per region | Fix to one row per object |

---

## Open questions

1. **Default object method.** `watershed` as the default for every stain except Abeta
   is pending a runtime comparison with `components` on one subject at level 0.
2. **Group maps.** Are per-region group means the only maps needed, or also voxel-wise
   group averages of the template-space field-fraction and count maps?
3. **Template-space view copy resolution.** Which `--template_spim_res` the explorer
   uses (one resolution, or a pyramid).
4. **MRI regional metrics and connectivity types that go in first** (also open in the
   portal plan).
