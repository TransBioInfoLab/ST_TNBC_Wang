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

  xy <- as.matrix(md[, c("x_um", "y_um"), drop = FALSE])
  median_nn_um <- if (nrow(xy) >= 2) {
    nn <- RANN::nn2(xy, xy, k = 2)
    median(nn$nn.dists[, 2], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  raw_counts <- get_assay_matrix(obj, assay = "RNA", layer = "counts")
  n_genes_detected <- sum(Matrix::rowSums(raw_counts > 0) > 0)

  data.frame(
    run_id = as.character(unique(md$run_id)[1]),
    count_id = as.character(unique(md$count_id)[1]),
    subarray_id = as.character(unique(md$subarray_id)[1]),
    subject_id = as.character(unique(md$subject_id)[1]),
    specimen_type = as.character(unique(md$specimen_type)[1]),
    baseline_primary_specimen = all(md$baseline_primary_specimen),
    counts_branch = as.character(unique(md$counts_branch)[1]),
    module_score_version = as.character(unique(md$module_score_version)[1]),
    
    n_spots = nrow(md),
    n_genes = nrow(raw_counts),
    n_genes_detected = n_genes_detected,
    n_annotation_supported = sum(md$manual_annotation_supported, na.rm = TRUE),

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
    "^subject_id$",
    "^specimen_type$",
    "^baseline_primary_specimen$",
    "^run_id$",
    "^counts_branch$",
    "^module_score_version$",
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

append_status <- function(status_row, status_file, run_id) {
  status_row <- data.frame(
    run_id = run_id,
    status_row,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
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

  details_file <- file.path(config$output_dir, "module_coverage",
    paste0(subarray_id, "_scoring_details.rds")
  )
  
  if (any(file.exists(c(required_existing, details_file)))) {
    stop(
      "Existing outputs found for ", subarray_id,
      "; use a fresh output directory."
    )
  }

  safe_message("Processing subarray: ", subarray_id)

  result <- extract_wang_subarray(
    parent = parent,
    subarray_id = subarray_id,
    min_annotation_pixels = config$min_annotation_biological_pixels
  )
  
  obj <- result$seurat
  obj <- perform_wang_qc(obj, min_umi = config$min_umi,
    max_artifact_fraction = config$max_artifact_fraction
  )
  obj <- add_wang_micron_coordinates(obj, spot_pitch_um = config$spot_pitch_um)
  obj <- normalize_wang_subarray(obj)
  obj <- compute_module_scores(obj, gene_sets,assay = "RNA", layer = "data")
  
  obj$run_id <- config$run_id
  obj$counts_branch <- config$counts_branch
  obj$module_score_version <- obj@misc$wang_module_scoring$version
  
  coverage <- calculate_module_coverage(obj, gene_sets, subarray_id, count_id)
  coverage$run_id <- config$run_id
  
  qc <- summarize_wang_qc(obj)
  spot_table <- make_spot_table(obj)
  
  normalized_expression <- get_assay_matrix(obj, assay = "RNA", layer = "data")
  
  stopifnot(
    identical(
      colnames(normalized_expression),
      spot_table$spot_id
    ),
    !anyDuplicated(spot_table$spot_id)
  )
  
  # Preserve mapping details even if Seurat objects are not saved.
  saveRDS(
    list(
      run_id = config$run_id,
      feature_map = obj@misc$wang_feature_map,
      module_scoring = obj@misc$wang_module_scoring,
      source_files = list(
        counts = normalizePath(
          file.path(config$counts_dir, paste0("TNBC", count_id, ".RDS")),
          mustWork = TRUE
        ),
        images = normalizePath(
          file.path(config$images_dir, paste0("TNBC", count_id, ".RDS")),
          mustWork = TRUE
        )
      )
    ),
    details_file,
    compress = TRUE
  )
  
  data.table::fwrite(coverage, coverage_file)
  data.table::fwrite(qc, qc_file)
  
  saveRDS(spot_table, spot_rds_file, compress = TRUE)
  
  if (isTRUE(config$save_spot_csv)) {
    data.table::fwrite(spot_table, spot_csv_file)
  }
  
  saveRDS(normalized_expression, normalized_file, compress = TRUE)
  
  if (isTRUE(config$save_seurat_rds)) {
    saveRDS(obj, object_file, compress = TRUE)
  }
  
  n_spots <- ncol(obj)
  
  rm(
    obj, result, spot_table,
    normalized_expression, coverage, qc
  )
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
  used_config_file <- normalizePath(config_file, mustWork = TRUE)
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(config$project_dir)
  
  config$output_dir <- config$working_output_dir
  
  if (file.exists(config$output_dir) &&
      !dir.exists(config$output_dir)) {
    stop("Output path is an existing file.")
  }
  
  if (dir.exists(config$output_dir) &&
      length(list.files(
        config$output_dir,
        all.files = TRUE,
        no.. = TRUE
      ))) {
    stop(
      "Output directory must be empty: ",
      config$output_dir
    )
  }
  
  threshold <- config$min_annotation_biological_pixels
  
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      !is.finite(threshold) ||
      threshold < 0) {
    stop(
      "Set min_annotation_biological_pixels ",
      "to one nonnegative number."
    )
  }
  
  config$counts_branch <- basename(
    normalizePath(config$counts_dir, mustWork = TRUE)
  )
  
  if (config$counts_branch != "countsNonCorrected") {
    stop(
      "This runner expects the ",
      "countsNonCorrected raw-count branch."
    )
  }
  
  if (is.null(config$run_id)) {
    config$run_id <- paste0(
      "wang_", config$release_version, "_",
      format(
        Sys.time(), "%Y%m%dT%H%M%S",
        tz = "UTC"
      ),
      "Z"
    )
  }
  
  if (!is.character(config$run_id) ||
      length(config$run_id) != 1L ||
      is.na(config$run_id) ||
      !nzchar(config$run_id)) {
    stop("run_id must be one nonempty character value.")
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

  if (!length(parent_ids)) {stop("No specimens selected.")}
  
  if (anyDuplicated(rownames(ids)) || anyNA(ids$id)) {
    stop("Invalid specimen/subarray identifiers.")
  }
  
  annotation_ids <- unique(as.character(
    ids$id[ids$hasAnnot %in% TRUE & as.character(ids$id) %in% parent_ids]
  ))
  
  if (length(annotation_ids) && is.null(config$annotation_dir)) {
    stop("annotation_dir is required for annotated specimens.")
  }
  
  code_files <- c(
    config = used_config_file,
    runner = file.path(config$project_dir, "02_run_wang_full_cohort.R"),
    helpers = config$function_file,
    gene_sets = config$gene_sets_file
  )
  
  annotation_files <- if (length(annotation_ids)) {
    setNames(
      file.path(
        config$annotation_dir,
        paste0("TNBC", annotation_ids, ".RDS")
      ),
      paste0("annotation_", annotation_ids)
    )
  } else {
    character()
  }
  
  input_files <- c(
    code_files,
    clinical = config$clinical_file,
    ids = config$ids_file,
    classification = config$classification_file,
    setNames(
      file.path(
        config$counts_dir,
        paste0("TNBC", parent_ids, ".RDS")
      ),
      paste0("counts_", parent_ids)
    ),
    setNames(
      file.path(
        config$images_dir,
        paste0("TNBC", parent_ids, ".RDS")
      ),
      paste0("images_", parent_ids)
    ),
    annotation_files
  )
  
  if (any(!file.exists(input_files))) {
    stop(
      "Missing inputs: ",
      paste(input_files[!file.exists(input_files)], collapse = ", ")
    )
  }
  
  input_manifest <- data.frame(
    run_id = config$run_id,
    input = names(input_files),
    path = normalizePath(unname(input_files), mustWork = TRUE),
    md5 = unname(tools::md5sum(input_files)),
    stringsAsFactors = FALSE
  )
  
  if (anyNA(input_manifest$md5)) {
    stop("Could not checksum every input.")
  }
  
  make_output_dirs(config$output_dir)
  
  logs_dir <- file.path(config$output_dir, "logs")
  code_dir <- file.path(logs_dir, "source_code")
  dir.create(code_dir, showWarnings = FALSE)
  
  if (!all(file.copy(
    code_files,
    file.path(code_dir, paste0(names(code_files), ".R")),
    overwrite = FALSE
  ))) {
    stop("Could not save all source-code snapshots.")
  }
  
  config$git_commit <- NA_character_
  git <- Sys.which("git")
  
  if (nzchar(git)) {
    commit <- tryCatch(
      suppressWarnings(system2(
        git,
        c("-C", shQuote(config$project_dir), "rev-parse", "HEAD"),
        stdout = TRUE,
        stderr = FALSE
      )),
      error = function(e) character()
    )
    
    if (length(commit) == 1L && is.null(attr(commit, "status"))) {
      config$git_commit <- commit[[1L]]
    }
  }
  
  saveRDS(config, file.path(logs_dir, "processing_config.rds"))
  saveRDS(
    list(gene_sets = GENE_SETS, alias_groups = wang_gene_alias_groups()),
    file.path(logs_dir, "scoring_definitions.rds")
  )
  
  data.table::fwrite(
    input_manifest,
    file.path(logs_dir, "processing_input_manifest.csv")
  )
  
  capture.output(
    sessionInfo(),
    file = file.path(
      logs_dir, "processing_session_info.txt"
    )
  )
  
  status_file <- file.path(
    config$output_dir,
    "status",
    "subarray_processing_status.csv"
  )

  for (count_id in parent_ids) {
    safe_message("Loading parent count ID: ", count_id)

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
          append_status(one_status, status_file, run_id = config$run_id)
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
            status_file,
            run_id = config$run_id
          )
        }
      }
    )
  }

  status_all <- data.table::fread(status_file)
  
  if (!all(c("run_id", "success") %in% names(status_all)) ||
      anyNA(status_all$success) ||
      !all(status_all$run_id %in% config$run_id)) {
    stop("Status records contain a different run_id.")
  }
  
  status_latest <- status_all[, .SD[.N], by = subarray_id]
  expected_ids <- rownames(ids)[as.character(ids$id) %in% parent_ids]
  
  if (!setequal(status_latest$subarray_id, expected_ids)) {
    stop(
      "Status records do not cover exactly ",
      "the selected subarrays."
    )
  }
  
  data.table::fwrite(
    status_latest,
    file.path(
      config$output_dir,
      "status",
      "subarray_processing_status_latest.csv"
    )
  )
  
  good_ids <- status_latest[success %in% TRUE, subarray_id]
  if (!length(good_ids)) {
    stop("No subarrays completed successfully.")
  }
  
  read_current_summary <- function(path, sid) {
    x <- data.table::fread(path)
    
    if (!all(c("run_id", "subarray_id") %in% names(x)) ||
    !nrow(x) ||
    !all(x$run_id %in% config$run_id) ||
    !identical(unique(as.character(x$subarray_id)), sid)) {
      stop(
        "Summary does not belong to the ",
        "current array/run: ", path
      )
    }
    x
  }
  
  qc_all <- data.table::rbindlist(Map(
    read_current_summary,
    file.path(config$output_dir, "qc_summaries", paste0(good_ids, "_qc.csv")),
    good_ids
  ))
  
  coverage_all <- data.table::rbindlist(Map(
    read_current_summary,
    file.path(config$output_dir, "module_coverage", paste0(good_ids, "_module_coverage.csv")),
    good_ids
  ))
  
  if (nrow(qc_all) != length(good_ids) ||
      anyDuplicated(qc_all$subarray_id) ||
      nrow(coverage_all) !=
      length(good_ids) * length(GENE_SETS) ||
      anyDuplicated(
        coverage_all[, .(subarray_id, module)]
      )) {
    stop("Duplicate or incomplete QC/module-summary records.")
  }
  
  data.table::setorder(qc_all, subarray_id)
  data.table::setorder(coverage_all, subarray_id, module)
  
  data.table::fwrite(
    qc_all,
    file.path(config$output_dir, "cohort_qc_summary.csv")
  )
  
  data.table::fwrite(
    coverage_all,
    file.path(
      config$output_dir,
      "cohort_module_coverage.csv"
    )
  )
  
  # Accept only the already-audited empty array and matching error.
  expected_empty <- status_latest[
    success %in% FALSE &
      subarray_id == "CN17_E2" &
      !is.na(error) &
      grepl(
        "No spots found for subarray: CN17_E2",
        error,
        fixed = TRUE
      ),
    subarray_id
  ]
  
  unexpected <- status_latest[success %in% FALSE &!subarray_id %in% expected_empty]
  
  full_run <- setequal(parent_ids, unique(as.character(ids$id)))
  
  accounting_ok <- !full_run || (
    nrow(status_latest) == 281L &&
      nrow(qc_all) == 280L &&
      sum(qc_all$n_spots) == 270310L &&
      data.table::uniqueN(qc_all$subject_id) == 92L
  )
  
  run_complete <- nrow(unexpected) == 0L && accounting_ok
  
  data.table::fwrite(
    data.frame(
      run_id = config$run_id,
      run_status = if (run_complete) {"complete"} else {"failed"},
      full_cohort = full_run,
      n_successful = length(good_ids),
      n_expected_empty = length(expected_empty),
      n_unexpected_failures = nrow(unexpected),
      n_spots = sum(qc_all$n_spots),
      accounting_ok = accounting_ok
    ),
    file.path(
      config$output_dir,
      "logs",
      "processing_run_summary.csv"
    )
  )
  
  if (!run_complete) {
    stop(
      "Run failed acceptance checks; inspect ",
      "the status and run-summary files."
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
