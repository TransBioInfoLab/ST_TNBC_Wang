#!/usr/bin/env Rscript

################################################################################
# 04_prepare_wang_release.R
#
# This script does not recreate normalized expression, spot metadata, or
# Seurat objects. Those files are created directly by Step 02.
#
# It adds release documentation to outputs_wang:
#   - cohort manifest
#   - file manifest with checksums
#   - data dictionary
#   - module-gene table
#   - processed-data README
#   - session information
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

find_wang_config <- function() {
  env_file <- Sys.getenv("WANG_TNBC_CONFIG", unset = "")
  if (nzchar(env_file)) return(env_file)

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_file <- sub("^--file=", "", file_arg[1])
    return(file.path(dirname(normalizePath(script_file)), "00_wang_config.R"))
  }

  file.path(getwd(), "00_wang_config.R")
}

config_file <- find_wang_config()

if (!file.exists(config_file)) {
  stop(
    "Configuration file not found: ", config_file,
    "\nRun the script from the TNBC_Wang folder or set WANG_TNBC_CONFIG."
  )
}

source(config_file)


describe_variable <- function(variable) {
  category <- dplyr::case_when(
    grepl("^clinical_", variable) ~ "clinical",
    grepl("^hist_score_", variable) ~ "histomorphology_score",
    grepl("^manual_fraction_", variable) ~ "manual_annotation_fraction",
    grepl("^manual_pixels_", variable) ~ "manual_annotation_pixels",
    grepl("^manual_", variable) ~ "manual_annotation_qc",
    grepl("^score_", variable) ~ "marker_module_score",
    grepl("^qc_", variable) ~ "qc_flag",
    grepl("^percent\\.", variable) ~ "qc_metric",

    variable %in% c(
      "x", "y", "new_x", "new_y",
      "pixel_x", "pixel_y", "x_um", "y_um"
      ) ~ "coordinate",

    variable %in% c(
      "spot_id", "study", "platform", "count_id",
      "subarray_id", "parent_slide", "subslide",
      "sample_label", "has_manual_annotation"
    ) ~ "identifier",

    TRUE ~ "other"
  )

  units <- dplyr::case_when(
    variable %in% c("x_um", "y_um") ~ "microns",
    variable %in% c("pixel_x", "pixel_y") ~ "image pixels",
    grepl("^manual_fraction_", variable) ~ "fraction from 0 to 1",
    grepl("^manual_pixels_", variable) ~ "pixels",
    grepl("^percent\\.", variable) ~ "percent",
    grepl("^score_", variable) ~ "expression-enrichment score",
    grepl("^hist_score_", variable) ~ "model score",
    TRUE ~ NA_character_
  )

  availability <- dplyr::case_when(
    grepl("^manual_", variable) ~ "annotated subarrays only or NA otherwise",
    TRUE ~ "processed subarrays"
  )

  description <- dplyr::case_when(
    variable == "spot_id" ~
      "Unique spot identifier; matches normalized-expression matrix column names.",

    variable == "count_id" ~
      "Parent source-sample/count-object identifier.",

    variable == "subarray_id" ~
      "Unique legacy-ST subarray identifier.",

    variable == "x_um" ~
      "Centered x coordinate in microns, derived from the original array coordinates.",

    variable == "y_um" ~
      "Centered y coordinate in microns, derived from the original array coordinates.",

    variable == "pixel_x" ~
      "Image x coordinate supplied by the source repository.",

    variable == "pixel_y" ~
      "Image y coordinate supplied by the source repository.",

    variable == "nCount_RNA" ~
      "Total UMI count for the spot.",

    variable == "nFeature_RNA" ~
      "Number of detected genes for the spot.",

    variable == "has_manual_annotation" ~
      "Whether the subarray has manual histopathology annotations.",

    variable == "qc_pass" ~
      "Combined UMI and artifact QC flag. The source-retained spots are not removed by this pipeline.",

    grepl(
      "^hist_score_",
      variable
    ) ~
      "Continuous source histomorphology-classification score; not an exact tissue fraction.",

    grepl(
      "^manual_fraction_",
      variable
    ) ~
      "Fraction of biologically annotated pixels assigned to the named manual category.",

    grepl(
      "^manual_pixels_",
      variable
    ) ~
      "Raw number of annotation pixels assigned to the named category within the spot.",

    grepl(
      "^score_.*_z$",
      variable
    ) ~
      "Within-subarray robust z-score of the corresponding marker-module mean.",

    grepl(
      "^score_.*_pos$",
      variable
    ) ~
      "Within-subarray robust z-score truncated at zero; not a cell proportion.",

    grepl("^score_", variable) ~
      "Mean log-normalized expression of the named marker module; a marker-module enrichment score, not a deconvolved cell-type proportion.",

    grepl(
      "^clinical_.*_(time|status)$",
      variable
    ) ~
      "Scalar component expanded from a source survival outcome.",

    grepl("^clinical_", variable) ~
      "Source sample-level clinical variable repeated across spots.",

    TRUE ~
      "See source code and source-repository documentation."
  )

  data.frame(
    variable = variable,
    category = category,
    units_or_scale = units,
    availability = availability,
    description = description,
    stringsAsFactors = FALSE
  )
}

build_metadata_dictionary <- function(spot_files) {
  type_rows <- lapply(
    spot_files,
    function(file) {
      table <- readRDS(file)

      data.frame(
        variable = names(table),
        data_type = vapply(
          table,
          function(x) {paste(class(x), collapse = "/")},
          character(1)
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  type_table <-
    data.table::rbindlist(type_rows,fill = TRUE) |>
    as.data.frame() |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      data_type = paste(sort(unique(data_type)), collapse = ";"),
      .groups = "drop"
    )

  dictionary <-
    data.table::rbindlist(
      lapply(type_table$variable, describe_variable)) |>
    as.data.frame() |>
    dplyr::left_join(type_table, by = "variable") |>
    dplyr::select(
      variable, category, data_type, units_or_scale, availability, description
    ) |>
    dplyr::arrange(category, variable)

  dictionary
}

run_prepare_release <- function(config = CONFIG) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  setwd(config$project_dir)

  output_dir <- config$working_output_dir
  manifest_dir <- file.path(output_dir, "manifests")
  documentation_dir <- file.path(output_dir, "documentation")
  logs_dir <- file.path(output_dir, "logs")

  dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(documentation_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  source(config$gene_sets_file)

  ids <- readRDS(config$ids_file)
  ids$subarray_id <- rownames(ids)

  status_file <- file.path(output_dir, "status",
    "subarray_processing_status_latest.csv")

  if (!file.exists(status_file)) {
    stop(
      "Latest status file not found: ",
      status_file
    )
  }

  status <- data.table::fread(status_file) |>
    as.data.frame()
  qc <- data.table::fread(file.path(output_dir, "cohort_qc_summary.csv")) |>
    as.data.frame() 
  eligibility_file <- file.path(config$validation_results_dir,
    "subarray_eligibility_for_review.csv")

  eligibility <- if (file.exists(eligibility_file)) {
    data.table::fread(eligibility_file) |>
      as.data.frame()
  } else {
    data.frame(subarray_id = character())
  }
  
  # Keep the manifest concise. The eligibility table already contains the
  # QC table and module-coverage columns, so joining it in full would duplicate
  # fields and create `.x`/`.y` names. Select only the manifest-level flags.
  
  make_output_relative <- function(path_value) {
    path_value <- as.character(path_value)
    path_value <- gsub("\\\\", "/", path_value)
    
    prefix <- paste0(gsub("\\\\", "/", output_dir), "/")
    
    has_prefix <- !is.na(path_value) & startsWith(path_value, prefix)
    
    path_value[has_prefix] <- substring(path_value[has_prefix], nchar(prefix) + 1)
    path_value
  }
  
  status_manifest <- status |>
    dplyr::transmute(
      count_id = as.character(count_id),
      subarray_id = as.character(subarray_id),
      processing_timestamp = as.character(timestamp),
      processing_success = as.logical(success),
      processing_skipped = as.logical(skipped),
      
      normalized_expression_file = make_output_relative(normalized_expression_file),
      spot_table_rds_file = make_output_relative(spot_table_rds_file),
      spot_table_csv_file = make_output_relative(spot_table_csv_file),
      processed_object_file = make_output_relative(object_file),
      processing_error = as.character(error)
    )
  
  qc_manifest <- qc |>
    dplyr::transmute(
      count_id = as.character(count_id),
      subarray_id = as.character(subarray_id),
      n_spots = as.integer(n_spots),
      n_genes = as.integer(n_genes),
      
      median_nCount_RNA,
      q05_nCount_RNA,
      q95_nCount_RNA,
      
      median_nFeature_RNA,
      q05_nFeature_RNA,
      q95_nFeature_RNA,
      
      median_percent_mt,
      q95_percent_mt,
      median_percent_ribo,
      median_percent_mtrnr,
      
      pct_umi_ge_500,
      pct_artifact_pass,
      pct_qc_pass,
      pct_classification_complete,
      
      median_nn_um
    )
  
  eligibility_manifest <- eligibility |>
    dplyr::select(
      dplyr::any_of(
        c("subarray_id","flag_no_spots","flag_lt100",
          "flag_100_199","eligible_primary","eligible_tumor",
          "eligible_lymphocyte","eligible_tls","review_status"
        )
      )
    ) |>
    dplyr::mutate(
      subarray_id = as.character(subarray_id)
    )
  

  cohort_manifest <- ids |>
    dplyr::transmute(
      count_id = as.character(id),
      subarray_id = as.character(subarray_id),
      parent_slide = as.character(slide),
      subslide = as.character(subslide),
      sample_label = as.character(ID),
      source_has_manual_annotation = as.logical(hasAnnot)
    ) |>
    dplyr::left_join(
      status_manifest,
      by = c("count_id","subarray_id")
    ) |>
    dplyr::left_join(
      qc_manifest,
      by = c("count_id", "subarray_id"),
      suffix = c("", "_qc")
    ) |>
    dplyr::left_join(eligibility_manifest, by = "subarray_id") |>
    dplyr::mutate(
      has_retained_spots = processing_success %in% TRUE,

      inclusion_status =
        dplyr::case_when(
          processing_success %in% TRUE ~ "included",
          grepl("No spots found", processing_error, fixed = TRUE) ~ "excluded_no_retained_spots",
          processing_success %in% FALSE ~ "processing_failed",
          TRUE ~ "not_processed"
        ),

      exclusion_or_qc_note =
        dplyr::case_when(
          inclusion_status == "excluded_no_retained_spots" ~
            "No retained spots in the source filtered count object.",
          subarray_id == "CN26_D1" ~
            "Included in the processed dataset; only 8 retained spots and not recommended for primary spatial analysis.",
          !is.na(review_status) & review_status != "Primary eligible" ~ review_status,
          inclusion_status == "processing_failed" ~ processing_error,
          TRUE ~  NA_character_
        ),

      release_version = config$release_version,
      release_date = config$release_date
    ) |>
    dplyr::arrange(suppressWarnings(as.numeric(count_id)), subarray_id)

  data.table::fwrite(
    cohort_manifest,
    file.path(manifest_dir, "cohort_manifest.csv")
  )

  module_gene_sets <-
    data.table::rbindlist(
      lapply(
        names(GENE_SETS),
        function(module_name) {
          data.frame(
            module = module_name,
            gene = GENE_SETS[[module_name]],
            interpretation =
              "Marker-module expression enrichment; not a deconvolved cell-type proportion.",
            stringsAsFactors = FALSE
          )
        }
      )
    )

  data.table::fwrite(
    module_gene_sets,
    file.path(documentation_dir, "module_gene_sets.csv")
  )

  spot_files <- list.files(
    file.path(output_dir, "spot_tables"),
    pattern = "_spots\\.rds$",
    full.names = TRUE
  )

  if (length(spot_files) == 0) {
    stop("No spot-table RDS files found.")
  }

  data_dictionary <- build_metadata_dictionary(spot_files)

  data.table::fwrite(
    data_dictionary,
    file.path(documentation_dir, "data_dictionary.csv")
  )

  processed_readme <- c(
    "# Wang TNBC ST processed dataset",
    "",
    paste0("Release version: ", config$release_version),
    "",
    "## Scope",
    "",
    "This release contains normalized expression, spot metadata, micron-scale coordinates, histomorphology scores, manual annotation fractions, marker-module enrichment scores, QC summaries, module-coverage summaries, and marker-score validation summaries.",
    "",
    "It does not contain reference-based cell-type deconvolution or finalized tumor-immune proximity outcomes.",
    "",
    "T-cell, CD8/cytotoxic, B-cell, plasma-cell, macrophage, and M2-macrophage variables are marker-module enrichment scores, not deconvolved cell-type proportions.",
    "",
    "## Main folders",
    "",
    "- `normalized_expression/`: sparse genes-by-spots log-normalized matrices.",
    "- `spot_tables/`: spot-level metadata in RDS and CSV.GZ formats.",
    "- `processed_objects/`: Seurat objects.",
    "- `qc_summaries/`: one QC summary per subarray.",
    "- `module_coverage/`: one module-coverage summary per subarray.",
    "- `abundance_validation/`: marker-score validation against manual annotations.",
    "- `manifests/`: cohort and file manifests.",
    "- `documentation/`: data dictionary and module definitions.",
    "",
    "## Normalization",
    "",
    "`log(1 + 10000 * raw_count / total_raw_counts_for_the_spot)`",
    "",
    "Rows are genes and columns are spot IDs. Matrix column names exactly match `spot_id` in the corresponding spot table.",
    "",
    "## Coordinates",
    "",
    "`x_um` and `y_um` are centered micron-scale coordinates based on the original regular array coordinates and a 150-micrometer center-to-center spot pitch.",
    "",
    "## Cohort accounting",
    "",
    "The source manifest contains 281 listed subarrays. The expected processed release contains 280 subarrays and 270,310 retained spots. CN17_E2 has no retained spots in the source filtered count object.",
    "",
    "## Important interpretation",
    "",
    "`hist_score_*` variables are continuous histomorphology-classification scores, not exact tissue fractions.",
    "",
    "`manual_fraction_*` variables are manual annotation fractions among biologically annotated pixels and are available only in annotated subarrays.",
    "",
    "`score_*` variables are marker-module expression-enrichment scores. They are not cell proportions."
  )

  writeLines(
    processed_readme,
    file.path(documentation_dir, "README_processed_data.md")
  )

  capture.output(
    sessionInfo(),
    file = file.path(logs_dir, "session_info.txt")
  )

  # Build the file manifest last so that it includes all documentation.
  existing_manifest <- file.path(manifest_dir, "file_manifest.csv")

  if (file.exists(existing_manifest)) {
    file.remove(existing_manifest)
  }

  files <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )

  file_info <- file.info(files)

  root_normalized <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  files_normalized <- normalizePath(files, winslash = "/", mustWork = TRUE)

  relative_paths <- substring(
    files_normalized,
    nchar(root_normalized) + 2
  )

  file_manifest <- data.frame(
    relative_path = relative_paths,
    file_name = basename(files),
    file_extension = tools::file_ext(files),
    file_size_bytes = file_info$size,
    md5 = unname(tools::md5sum(files)),
    release_version = config$release_version,
    release_date = config$release_date,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(file_manifest, existing_manifest)

  message(
    "Release manifests and documentation created in: ",
    output_dir
  )

  invisible(
    list(
      cohort_manifest = cohort_manifest,
      data_dictionary = data_dictionary,
      file_manifest = file_manifest
    )
  )
}

if (sys.nframe() == 0) {
  release_results <- run_prepare_release(CONFIG)
}
