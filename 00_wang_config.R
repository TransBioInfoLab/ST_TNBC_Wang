################################################################################
# 00_wang_config.R
#
# Local configuration for the Wang TNBC ST pipeline.
# All pipeline scripts are stored directly in the TNBC_Wang project folder.
################################################################################

CONFIG <- list(
  # Edit only this path if the project is moved.
  project_dir = "/Users/zhixinmao/Desktop/st_analysis/ST_TNBC_Wang",

  # Code files in the project root
  function_file = "01_wang_preprocess_functions.R",
  gene_sets_file = "01_wang_gene_sets.R",

  # Source data
  clinical_file = "Clinical/Clinical.RDS",
  ids_file = "Clinical/ids.RDS",
  classification_file = "classification/classifAll.RDS",
  counts_dir = "Robjects/countsNonCorrected",
  images_dir = "Robjects/images",
  annotation_dir = "Robjects/annotsBySpot",

  working_output_dir = "outputs_wang",

  # Preprocessing settings
  min_umi = 500,
  max_artifact_fraction = 0.01,
  spot_pitch_um = 150,
  save_seurat_rds = TRUE,
  save_spot_csv = TRUE,

  # Use 1 and 3 for a pilot, then change both to NULL for the full run.
  start_parent = NULL,
  stop_parent = NULL,

  # outputs_wang is new, so TRUE is appropriate for the first corrected run.
  overwrite = TRUE,

  # Marker-score validation
  annotated_spots_file =
    file.path("outputs_wang", "abundance_validation", "annotated_spots.rds"),

  validation_results_dir =
    file.path("outputs_wang", "abundance_validation"),

  primary_min_spots = 200,
  sensitivity_min_spots = 100,
  manual_positive_threshold = 0.25,
  min_complete_spots_for_rho = 20,
  min_positive_for_auc = 5,
  min_negative_for_auc = 5,

  # Release/documentation version
  release_version = "1.0",
  release_date = format(Sys.Date(), "%Y-%m-%d")
)
