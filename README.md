# Wang TNBC spatial transcriptomics processing

All scripts are stored directly in the `ST_TNBC_Wang` project folder.

## Folder layout

```text
ST_TNBC_Wang/
├── 00_wang_config.R
├── 01_wang_gene_sets.R
├── 01_wang_preprocess_functions.R
├── 02_run_wang_full_cohort.R
├── 03_validate_wang_marker_scores.R
├── 04_prepare_wang_release.R
├── Clinical/
├── classification/
├── Robjects/
└── outputs_wang/
```

## First run: pilot

In `00_wang_config.R`, use:

```r
start_parent = 1
stop_parent = 3
overwrite = TRUE
```

Run:

```bash
Rscript 02_run_wang_full_cohort.R
```

Review the pilot outputs.

## Full corrected cohort

Delete the pilot `outputs_wang` folder or keep `overwrite = TRUE`.

Change:

```r
start_parent = NULL
stop_parent = NULL
```

Run:

```bash
Rscript 02_run_wang_full_cohort.R
```

The preprocessing step directly creates:

- normalized expression matrices;
- spot tables in RDS and CSV.GZ formats;
- Seurat objects;
- QC summaries;
- module-coverage summaries.

## Marker-score validation

Run:

```bash
Rscript 03_validate_wang_marker_scores.R
```

This validates tumor, broad lymphocyte, and TLS candidate scores against manual annotations.

Lineage-specific T-cell, CD8/cytotoxic, B-cell, plasma-cell, macrophage, and M2 variables remain marker-module enrichment scores, not deconvolved proportions.

## Prepare documentation

Run:

```bash
Rscript 04_prepare_wang_release.R
```

This does not recreate the biological outputs. It creates manifests, checksums, the data dictionary, module definitions, and release documentation inside `outputs_wang`.


## GitHub versus Dropbox

Private GitHub repository:

- configuration template;
- gene sets;
- all R scripts;
- this README.

Dropbox:

- upload the complete verified `outputs_wang` folder.
