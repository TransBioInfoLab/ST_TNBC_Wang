#!/usr/bin/env Rscript

################################################################################
# 02_run_wang_full_cohort.R
#
# Corrected raw-data-to-processed-data pipeline.
# Writes a new clean run to TNBC_Wang/outputs_wang.
# TNBC_Wang/outputs_wang_old is never modified.
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(RANN)
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


safe_message <- function(...) {
  message(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste0(..., collapse = "")
  )
}

make_output_dirs <- function(output_dir) {
  dirs <- file.path(
    output_dir,
    c(
      "processed_objects",
      "normalized_expression",
      "spot_tables",
      "qc_summaries",
      "module_coverage",
      "abundance_validation",
      "manifests",
      "documentation",
      "logs",
      "status"
    )
  )

  invisible(
    lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  )
}

summarize_wang_qc <- function(obj) {
  md <- obj@meta.data

  qfun <- function(x, p) {
    as.numeric(
      stats::quantile(x, p, na.rm = TRUE, names = FALSE)
    )
  }

  xy <- as.matrix(
    md[, c("x_um", "y_um"), drop = FALSE]
  )

  median_nn_um <- if (nrow(xy) >= 2) {
    nn <- RANN::nn2(xy, xy, k = 2)
    median(nn$nn.dists[, 2], na.rm = TRUE)
  } else {
    NA_real_
  }

  data.frame(
    count_id = as.character(unique(md$count_id)[1]),
    subarray_id = as.character(unique(md$subarray_id)[1]),

    n_spots = nrow(md),
    n_genes = nrow(obj),

    median_nCount_RNA = median(md$nCount_RNA, na.rm = TRUE),
    q05_nCount_RNA = qfun(md$nCount_RNA, 0.05),
    q95_nCount_RNA = qfun(md$nCount_RNA, 0.95),
    median_nFeature_RNA = median(md$nFeature_RNA, na.rm = TRUE),
    q05_nFeature_RNA = qfun(md$nFeature_RNA, 0.05),
    q95_nFeature_RNA = qfun(md$nFeature_RNA, 0.95),
    median_percent_mt = median(md$percent.mt, na.rm = TRUE),
    q95_percent_mt = qfun(md$percent.mt, 0.95),
    median_percent_ribo = median(md$percent.ribo, na.rm = TRUE),
    median_percent_mtrnr = median(md$percent.mtrnr, na.rm = TRUE),
    pct_umi_ge_500 = mean(md$qc_pass_umi, na.rm = TRUE) * 100,
    pct_artifact_pass = mean(md$qc_pass_artifact, na.rm = TRUE) * 100,
    pct_qc_pass = mean(md$qc_pass, na.rm = TRUE) * 100,

    pct_classification_complete =
      if ("hist_score_Tumor" %in% colnames(md)) {
        mean(is.finite(md$hist_score_Tumor)) * 100
      } else {
        0
      },

    has_manual_annotation =
      as.logical(unique(md$has_manual_annotation)[1]),

    median_nn_um = median_nn_um,

    stringsAsFactors = FALSE
  )
}

make_spot_table <- function(obj) {
  md <- obj@meta.data |>
    tibble::rownames_to_column("spot_id")

  keep_patterns <- c(
    "^spot_id$",
    "^study$",
    "^platform$",
    "^count_id$",
    "^subarray_id$",
    "^parent_slide$",
    "^subslide$",
    "^sample_label$",
    "^has_manual_annotation$",
    "^clinical_",
    "^x$",
    "^y$",
    "^new_x$",
    "^new_y$",
    "^pixel_x$",
    "^pixel_y$",
    "^x_um$",
    "^y_um$",
    "^coordinate_",
    "^microns_",
    "^assumed_",
    "^nCount_RNA$",
    "^nFeature_RNA$",
    "^percent\\.",
    "^qc_",
    "^hist_",
    "^manual_",
    "^score_"
  )

  keep_cols <- unique(
    unlist(
      lapply(
        keep_patterns,
        function(pattern) {
          grep(
            pattern,
            colnames(md),
            value = TRUE
          )
        }
      )
    )
  )

  out <- md[, keep_cols, drop = FALSE]

  validate_scalar_spot_table(
    out,
    table_name = paste0(
      "Spot table for ",
      unique(out$subarray_id)[1]
    )
  )

  out
}

append_status <- function(status_row, status_file) {
  data.table::fwrite(
    status_row,
    status_file,
    append = file.exists(status_file),
    col.names = !file.exists(status_file)
  )
}

process_wang_subarray <- function(parent, subarray_id, config, gene_sets) {
  count_id <- as.character(parent$count_id)

  spot_rds_rel <- file.path("spot_tables", paste0(subarray_id, "_spots.rds"))
  spot_csv_rel <- file.path("spot_tables", paste0(subarray_id, "_spots.csv.gz"))
  normalized_rel <- file.path("normalized_expression",
    paste0(subarray_id, "_lognormalized_expression.rds"))
  object_rel <- file.path("processed_objects", paste0(subarray_id, "_seurat.rds"))
  spot_rds_file <- file.path(config$output_dir, spot_rds_rel)
  spot_csv_file <- file.path(config$output_dir, spot_csv_rel)
  normalized_file <- file.path(config$output_dir, normalized_rel)
  object_file <- file.path(config$output_dir, object_rel)
  qc_file <- file.path(config$output_dir, "qc_summaries",
    paste0(subarray_id, "_qc.csv"))
  coverage_file <- file.path(config$output_dir, "module_coverage",
    paste0(subarray_id, "_module_coverage.csv"))

  required_existing <- c(spot_rds_file, normalized_file, qc_file, coverage_file)

  if (isTRUE(config$save_spot_csv)) {
    required_existing <- c(
      required_existing,
      spot_csv_file
    )
  }

  if (isTRUE(config$save_seurat_rds)) {
    required_existing <- c(
      required_existing,
      object_file
    )
  }

  complete_existing <- all(file.exists(required_existing))

  if (!isTRUE(config$overwrite) && complete_existing) {
    old_spots <- readRDS(spot_rds_file)

    safe_message("Skipping existing subarray: ", subarray_id)

    return(
      data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        count_id = count_id,
        subarray_id = subarray_id,
        success = TRUE,
        skipped = TRUE,
        n_spots = nrow(old_spots),

        normalized_expression_file = normalized_rel,
        spot_table_rds_file = spot_rds_rel,

        spot_table_csv_file =
          if (isTRUE(config$save_spot_csv)) {spot_csv_rel} else {NA_character_},

        object_file =
          if (isTRUE(config$save_seurat_rds)) {object_rel} else {NA_character_},

        error = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  safe_message(
    "Processing subarray: ",
    subarray_id
  )

  result <- extract_wang_subarray(parent = parent, subarray_id = subarray_id)

  obj <- result$seurat
  obj <- perform_wang_qc(obj, min_umi = config$min_umi,
    max_artifact_fraction = config$max_artifact_fraction)
  obj <- add_wang_micron_coordinates(obj, spot_pitch_um = config$spot_pitch_um)
  obj <- normalize_wang_subarray(obj)
  obj <- compute_module_scores(obj, gene_sets, assay = "RNA", layer = "data")
  coverage <- calculate_module_coverage(obj, gene_sets, subarray_id, count_id)
  data.table::fwrite(coverage, coverage_file)
  qc <- summarize_wang_qc(obj)

  data.table::fwrite(qc, qc_file)
  spot_table <- make_spot_table(obj)

  saveRDS(spot_table, spot_rds_file,compress = TRUE)

  if (isTRUE(config$save_spot_csv)) {
    data.table::fwrite(
      spot_table,
      spot_csv_file
    )
  }

  normalized_expression <- get_assay_matrix(obj, assay = "RNA",layer = "data")

  stopifnot(
    identical(
      colnames(normalized_expression),
      spot_table$spot_id
    )
  )

  saveRDS(normalized_expression, normalized_file, compress = TRUE)

  if (isTRUE(config$save_seurat_rds)) {
    saveRDS(obj, object_file, compress = TRUE)
  }

  n_spots <- ncol(obj)

  rm(obj, result, spot_table, normalized_expression, coverage, qc)

  gc(verbose = FALSE)

  data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    count_id = count_id,
    subarray_id = subarray_id,
    success = TRUE,
    skipped = FALSE,
    n_spots = n_spots,
    normalized_expression_file = normalized_rel,
    spot_table_rds_file = spot_rds_rel,
    spot_table_csv_file =
      if (isTRUE(config$save_spot_csv)) {
        spot_csv_rel
      } else {
        NA_character_
      },
    object_file =
      if (isTRUE(config$save_seurat_rds)) {
        object_rel
      } else {
        NA_character_
      },
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_wang_full_cohort <- function(config = CONFIG) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  setwd(config$project_dir)

  config$output_dir <- config$working_output_dir

  make_output_dirs(config$output_dir)

  if (!file.exists(config$function_file)) {
    stop(
      "Function file not found: ",
      config$function_file
    )
  }

  if (!file.exists(config$gene_sets_file)) {
    stop(
      "Gene-set file not found: ",
      config$gene_sets_file
    )
  }

  source(config$function_file)
  source(config$gene_sets_file)

  if (!exists("GENE_SETS")) {
    stop(
      "GENE_SETS is not defined after ",
      "sourcing the gene-set file."
    )
  }

  clinical <- readRDS(config$clinical_file)

  ids <- readRDS(config$ids_file)

  classification <- readRDS(config$classification_file)

  if (!"id" %in% colnames(ids)) {
    stop(
      "ids does not contain column `id`."
    )
  }

  parent_ids <- unique(as.character(ids$id))

  parent_ids <- parent_ids[!is.na(parent_ids) & nzchar(parent_ids)]

  numeric_ids <- suppressWarnings(as.numeric(parent_ids))

  parent_ids <- if (all(is.finite(numeric_ids))) {
    parent_ids[order(numeric_ids)]
  } else {
    sort(parent_ids)
  }

  if (!is.null(config$start_parent)) {
    parent_ids <- parent_ids[as.numeric(parent_ids) >= as.numeric(config$start_parent)]
  }

  if (!is.null(config$stop_parent)) {
    parent_ids <- parent_ids[as.numeric(parent_ids) <= as.numeric(config$stop_parent)]
  }

  status_file <- file.path(config$output_dir, "status", "subarray_processing_status.csv")

  if (isTRUE(config$overwrite) && file.exists(status_file)) {
    file.remove(status_file)
  }

  for (count_id in parent_ids) {
    safe_message(
      "Loading parent count ID: ", count_id
    )

    tryCatch(
      {parent <- load_wang_parent(
          count_id = count_id,
          ids = ids,
          clinical = clinical,
          counts_dir = config$counts_dir,
          images_dir = config$images_dir,
          classification = classification,
          annotation_dir = config$annotation_dir
        )

        subarray_ids <- rownames(parent$manifest_rows)

        for (subarray_id in subarray_ids) {
          one_status <- tryCatch(
            process_wang_subarray(
              parent,
              subarray_id,
              config,
              GENE_SETS
            ),
            error = function(e) {
              safe_message( "ERROR in ", subarray_id, ": ", conditionMessage(e))

              data.frame(
                timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                count_id = as.character(count_id),
                subarray_id = as.character(subarray_id),
                success = FALSE,
                skipped = FALSE,
                n_spots = NA_integer_,
                normalized_expression_file = NA_character_,
                spot_table_rds_file = NA_character_,
                spot_table_csv_file = NA_character_,
                object_file = NA_character_,
                error = conditionMessage(e),
                stringsAsFactors = FALSE
              )
            }
          )

          append_status(one_status, status_file)
        }

        rm(parent)
        gc(verbose = FALSE)
      },
      error = function(e) {
        safe_message("PARENT-LEVEL ERROR for ",count_id,": ",conditionMessage(e))

        parent_rows <- ids[as.character(ids$id) == as.character(count_id), ,drop = FALSE]

        parent_subarrays <- rownames(parent_rows)

        if (length(parent_subarrays) == 0) {
          parent_subarrays <- NA_character_
        }

        for (
          subarray_id in
          parent_subarrays
        ) {
          append_status(
            data.frame(
              timestamp = format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
              count_id = as.character(count_id),
              subarray_id = as.character(subarray_id),
              success = FALSE,
              skipped = FALSE,
              n_spots = NA_integer_,
              normalized_expression_file = NA_character_,
              spot_table_rds_file = NA_character_,
              spot_table_csv_file = NA_character_,
              object_file = NA_character_,
              error = paste0("Parent-level error: ",conditionMessage(e)),
              stringsAsFactors = FALSE
            ),
            status_file
          )
        }
      }
    )
  }

  status_all <- data.table::fread(status_file)

  status_latest <- status_all[, .SD[.N], by = subarray_id]

  data.table::fwrite(
    status_latest,
    file.path(config$output_dir, "status",
      "subarray_processing_status_latest.csv"
    )
  )

  qc_files <- list.files(
    file.path(config$output_dir, "qc_summaries"),
    pattern = "_qc\\.csv$",
    full.names = TRUE
  )

  if (length(qc_files) > 0) {
    qc_all <- data.table::rbindlist(
      lapply(qc_files, data.table::fread),
      fill = TRUE
    )

    data.table::setorder(qc_all, subarray_id)

    data.table::fwrite(
      qc_all,
      file.path(config$output_dir, "cohort_qc_summary.csv")
    )
  }

  coverage_files <- list.files(
    file.path(config$output_dir, "module_coverage"),
    pattern = "_module_coverage\\.csv$",
    full.names = TRUE
  )

  if (length(coverage_files) > 0) {
    coverage_all <-
      data.table::rbindlist(
        lapply(coverage_files, data.table::fread),
        fill = TRUE
      )

    data.table::setorder(coverage_all, subarray_id, module)

    data.table::fwrite(
      coverage_all,
      file.path(
        config$output_dir,
        "cohort_module_coverage.csv"
      )
    )
  }

  safe_message(
    "Finished. Successful: ",
    sum(status_latest$success, na.rm = TRUE),
    "; failed/excluded: ",
    sum(status_latest$success %in% FALSE, na.rm = TRUE)
  )

  invisible(status_latest)
}

if (sys.nframe() == 0) {
  run_wang_full_cohort(CONFIG)
}
