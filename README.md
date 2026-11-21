# Wang TNBC spatial transcriptomics preprocessing

Version covered by this workflow: **v1.1.0**.

This repository contains the R code for corrected preprocessing, marker-score validation, representative H&E/spot overlay QA, and final assembly of the Wang TNBC spatial transcriptomics dataset. Raw inputs and processed data are distributed separately from this source-code repository.

Source study: Wang et al., *Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer with potential clinical implications*, Nature Communications (2024). [Study and data availability](https://doi.org/10.1038/s41467-024-54145-w).

## Cohort and release scope

The verified full-cohort release contains:

| Item | Number |
| --- | ---: |
| Listed arrays/sections | 281 |
| Successfully processed arrays/sections | 280 |
| Retained spots | 270,310 |
| Source specimens | 94 |
| Independent patients | 92 |
| Arrays checked for source-coordinate consistency | 280 |
| Representative arrays with recorded technical image/spot reviews | 12 |

`CN17_E2` has no retained spots and is the single expected empty array. Relapse specimen 30 belongs to patient 53, and relapse specimen 66 belongs to patient 58. Specimens and repeated sections must not be counted as independent patients.

The release covers preprocessing, annotation agreement and technical image/spot QA. Full-resolution tumor-mask registration, pathology boundary approval, endpoint-specific robustness and clinical prediction are separate analyses.

## Scripts

Keep these seven R files together in the project root:

| File | Purpose |
| --- | --- |
| [00_wang_config.R](00_wang_config.R) | Paths, preprocessing settings, validation thresholds and version |
| [01_wang_gene_sets.R](01_wang_gene_sets.R) | Fixed marker panels |
| [01_wang_preprocess_functions.R](01_wang_preprocess_functions.R) | Import, gene mapping, normalization, scoring and metadata helpers |
| [02_run_wang_full_cohort.R](02_run_wang_full_cohort.R) | Process source counts and create frozen provenance records |
| [03_validate_wang_marker_scores.R](03_validate_wang_marker_scores.R) | Validate scores against supported manual spot annotations |
| [04_check_wang_spatial_overlays.R](04_check_wang_spatial_overlays.R) | Check coordinates and generate representative H&E/spot overlays before assembly |
| [05_prepare_wang_release.R](05_prepare_wang_release.R) | Assemble the reviewed release and verify its full file inventory/checksums |

The old `04_prepare_wang_release.R` has been replaced by the current Steps 04 and 05. Remove that obsolete script when updating the repository so the step order is unambiguous. Files beginning with `00` and `01` are configuration/helpers sourced by the execution scripts.

## Requirements and local inputs

The audited environment used R 4.5.1. Required packages include `Seurat`, `SeuratObject`, `Matrix`, `data.table`, `dplyr`, `tibble`, `RANN`, `tidyr`, `ggplot2`, `pROC` and the Bioconductor package `EBImage`. Processing, validation, overlay and release session information is archived with the outputs.

Prepare the source objects locally using the following layout:

```text
ST_TNBC_Wang/
├── README.md
├── 00_wang_config.R
├── 01_wang_gene_sets.R
├── 01_wang_preprocess_functions.R
├── 02_run_wang_full_cohort.R
├── 03_validate_wang_marker_scores.R
├── 04_check_wang_spatial_overlays.R
├── 05_prepare_wang_release.R
├── Clinical/
│   ├── Clinical.RDS
│   └── ids.RDS
├── classification/
│   ├── classifAll.RDS
│   └── loosReg.RData
└── Robjects/
    ├── countsNonCorrected/TNBC*.RDS
    ├── images/TNBC*.RDS
    └── annotsBySpot/TNBC*.RDS
```

Edit `CONFIG$project_dir` in `00_wang_config.R` to point to this local project. Check all input paths. The source `countsNonCorrected` branch is required for the corrected preprocessing workflow.

The default processing directory is `outputs_wang_v1.1`. The validation paths must point into the same processing directory. Step 02 requires a **new, empty output directory**, even when `overwrite = TRUE`; preserve an earlier run and choose another directory for a new run. For the full cohort, set both `start_parent` and `stop_parent` to `NULL`. Partial processing pilots do not satisfy the full-cohort requirements for final assembly.

The configuration value `release_version = "1.1"` is canonicalized to final release version `1.1.0`.

## Run the workflow

Run commands from the project root.

### Step 02: full-cohort preprocessing

```bash
Rscript 02_run_wang_full_cohort.R
```

This creates sparse normalized expression, Seurat objects, RDS/CSV spot tables, QC summaries, module-coverage records and frozen code/configuration/input-checksum snapshots. Check that `logs/processing_run_summary.csv` reports `run_status = complete`, `full_cohort = TRUE` and zero unexpected failures.

### Step 03: marker-score validation

```bash
Rscript 03_validate_wang_marker_scores.R
```

Validation distinguishes full-data source predictions, original author specimen-held-out predictions and fixed marker-module annotation agreement. Main comparisons use supported spots and matched evaluable sections; available-predictor summaries also report coverage. Outputs are written to the configured `abundance_validation` directory.

### Step 04: representative spatial overlay QA

```bash
Rscript 04_check_wang_spatial_overlays.R
```

Step 04 reads completed processing outputs directly, without requiring a previously assembled release. It also recognizes the existing layout in which those outputs are inside `outputs_wang_v1.1/Wang_TNBC_v1.1.0`.

Each run creates a fresh directory under `spatial_overlay_QA_v1.1.0/run_<UTC timestamp>/`. Review the contact sheet and individual four-panel figures: original H&E, spot centers, supported manual Tumor fractions and UMI counts. Overlay panels use larger outlined markers and a lightened H&E background; these display settings do not change coordinates or expression values.

The default selection contains 12 arrays, including all four reflected source mappings, several rotation angles and annotation extremes. Complete `overlay_review.csv` using these status values:

| Column | Meaning and accepted decisions |
| --- | --- |
| `image_alignment_review` | `pass`, `fail`, `uncertain` or `pending`; use `pass` when the technical overlay check is satisfactory |
| `rotation_flip_review` | Same decisions; evaluate the displayed source transformation, without applying another rotation/reflection |
| `manual_annotation_review` | Same decisions, or `not_available` where manual annotation is unavailable; retain `pending` when not assessed |

Final assembly requires passing alignment/orientation decisions for every selected array, at least eight representative arrays, all reflected cases, and no unresolved automatic-coordinate or explicit manual-annotation failures. Pending manual-annotation assessment can remain documented separately from technical release readiness.

The editable review CSV is excluded from the Step 04 output checksum manifest. Generated figures and audit tables are protected; regenerate a QA run if they change. The final assembly independently hashes the reviewed CSV as part of the release.

### Step 05: final release assembly

An optional preflight check performs the numerical/provenance checks and reports readiness without assembling or replacing a release:

```bash
Rscript 05_prepare_wang_release.R --dry-run
```

To assemble the release:

```bash
Rscript 05_prepare_wang_release.R
```

By default, Step 05 selects the latest completed QA run matching the processing run. To select a specific reviewed run, supply its path:

```bash
Rscript 05_prepare_wang_release.R --overlay-dir=spatial_overlay_QA_v1.1.0/run_YYYYMMDDTHHMMSSZ
```

Replace the example run name with the actual directory produced by Step 04. The final folder remains **`outputs_wang_v1.1/Wang_TNBC_v1.1.0`**. Assembly occurs in a new staging directory. After successful verification, an existing final folder is retained in a sibling backup, and the staged folder takes its original name. No renormalization or rescoring occurs in Step 05.

Moved processing/validation inputs are resolved by relative paths and checked against their original hashes. The previously omitted frozen scoring-definitions file can be recovered from archived per-array details only when its original byte checksum is reproduced exactly; the file and recovery record are included in the new release.

## Final processed-data folder

```text
outputs_wang_v1.1/Wang_TNBC_v1.1.0/
├── normalized_expression/
├── processed_objects/
├── spot_tables/
├── qc_summaries/
├── module_coverage/
├── abundance_validation/
├── status/
├── logs/
├── documentation/
├── manifests/
├── verification/
│   └── spatial_overlays/
├── cohort_qc_summary.csv
└── cohort_module_coverage.csv
```

The selected QA run, archived overlay code, final review decisions and figures are included under `verification/spatial_overlays`. Source H&E objects are referenced by path/checksum rather than embedded in Seurat objects; the representative rendered figures are included.

## Measurement interpretation

- Expression uses `log(1 + 10000 * raw_count / spot_total)` from `countsNonCorrected`. Measured zero-expression genes are retained.
- Source-to-assay gene mapping and aliases are explicit. Marker scores use fixed panels and denominators; an unavailable required gene makes the corresponding score `NA` rather than silently shortening its panel.
- Marker-module scores measure expression enrichment, not cell proportions. Within-array `_z` and `_pos` scores describe local relative enrichment.
- Manual fractions use biological annotation pixels. Validation requires **strictly more than 1,000** such pixels; unsupported targets are not tumor-negative observations.
- `hist_score_*` predictions are expression-derived full-data fits. Their annotation agreement includes training data. `hist_loo_*` predictions use the authors' specimen-held-out scheme, which does not jointly hold out linked primary/relapse specimens from the same patient.
- The mean-expression `TLS_WANG_30` module does not reproduce every published weighted signature or establish external predictive performance.
- Spatial distances use the 150-micrometer lattice pitch and remain within arrays. Coordinates are not registered across sections. Gaps and partial capture require support-aware region selection; an overlay does not validate a precise tumor boundary or direct cell-to-cell contact.
- Clinical times/statuses remain as supplied. Endpoint definition, follow-up handling, patient-level resampling, raw-versus-author-corrected sensitivity and repeat-section agreement of the intended endpoint belong to subsequent analyses.
