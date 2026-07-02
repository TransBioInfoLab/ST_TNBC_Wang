#!/usr/bin/env Rscript

################################################################################
# 03_validate_wang_marker_scores.R
#
# Validates morphology and marker-module enrichment scores against manual
# annotation fractions in annotated subarrays.
#
# This script does not create deconvolved cell-type proportions.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pROC)
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


read_validation_spot_table <- function(file) {
  x <- readRDS(file)
  x <- as.data.frame(x, check.names = FALSE)

  fixed_cols <- c(
    "spot_id", "study", "platform", "count_id", "subarray_id",
    "parent_slide", "subslide", "sample_label", "x", "y", "new_x",
    "new_y", "pixel_x", "pixel_y", "x_um", "y_um", "nCount_RNA",
    "nFeature_RNA", "percent.mt", "percent.ribo", "percent.mtrnr",
    "qc_pass","has_manual_annotation"
  )

  pattern_cols <- grep(
    paste(
      c("^hist_score_",
        "^manual_fraction_",
        "^manual_pixels_",
        "^score_"),
      collapse = "|"
    ), names(x), value = TRUE)

  keep <- intersect(unique(c(fixed_cols, pattern_cols)), names(x))

  out <- x[, keep, drop = FALSE]

  scalar_ok <- vapply(
    out,
    function(z) {
      is.atomic(z) && is.null(dim(z)) && length(z) == nrow(out)},
    logical(1)
  )

  if (any(!scalar_ok)) {
    stop(
      "Non-scalar columns remain in ",
      basename(file), ": ",
      paste( names(out)[!scalar_ok], collapse = ", ")
    )
  }

  out
}

build_annotated_spots <- function(config = CONFIG) {
  output_dir <- config$working_output_dir
  qc_file <- file.path(output_dir, "cohort_qc_summary.csv")
  spot_dir <- file.path(output_dir, "spot_tables")
  qc <- data.table::fread(qc_file)

  annotated_ids <- qc |>
    dplyr::filter(has_manual_annotation) |>
    dplyr::pull(subarray_id)

  spot_files <- file.path(spot_dir, paste0(annotated_ids, "_spots.rds"))

  missing_files <- spot_files[
    !file.exists(spot_files)
  ]

  if (length(missing_files) > 0) {
    stop(
      "Missing annotated spot tables: ",
      paste(
        basename(missing_files),
        collapse = ", "
      )
    )
  }

  annotated_spots <- data.table::rbindlist(
      lapply(spot_files, read_validation_spot_table),
      fill = TRUE, use.names = TRUE
    ) |> as.data.frame()

  annotated_spots <- annotated_spots |>
    dplyr::mutate(
      manual_tumor = manual_fraction_Tumor + manual_fraction_Tumor.region,
      manual_lymphocyte = manual_fraction_Lymphocyte + 0.5 * manual_fraction_High.TIL.stroma,
      manual_tls = manual_fraction_Lymphoid.nodule,
      manual_tumor_positive = manual_tumor > config$manual_positive_threshold,
      manual_lymphocyte_positive = manual_lymphocyte > config$manual_positive_threshold,
      manual_tls_positive = manual_tls > config$manual_positive_threshold
    )

  dir.create(dirname(config$annotated_spots_file), recursive = TRUE, showWarnings = FALSE)

  saveRDS(annotated_spots, config$annotated_spots_file, compress = TRUE)

  annotated_spots
}

safe_cor <- function(x, y, min_n = 20) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < min_n) return(NA_real_)
  if (stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) {
    return(NA_real_)
  }
  
  suppressWarnings(
    stats::cor(x[keep], y[keep], method = "spearman")
  )
}

safe_auc <- function(x, y, min_pos = 5, min_neg = 5) {
  keep <- is.finite(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  
  if (sum(y == TRUE) < min_pos ||
    sum(y == FALSE) < min_neg ||
    length(unique(x)) < 2) {
    return(NA_real_)
  }
  
  ro <- tryCatch(
    pROC::roc(y, x, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(ro)) return(NA_real_)
  as.numeric(pROC::auc(ro))
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  
  as.numeric(
    stats::quantile(x, p, names = FALSE)
  )
}

rank01 <- function(x) {
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  
  if (!any(keep)) return(out)
  
  r <- rank(x[keep], ties.method = "average")
  
  out[keep] <- if (length(r) == 1) { 0.5 } else {
    (r - 1) / (length(r) - 1)
  }
  
  out
}

build_eligibility_table <- function(qc, coverage, config = CONFIG) {
  coverage_wide <- coverage |>
    dplyr::select(subarray_id, module, fraction_detected) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      names_from = module,
      values_from = fraction_detected,
      names_prefix = "coverage_"
    )
  
  qc |> dplyr::left_join(coverage_wide, by = "subarray_id") |>
    dplyr::mutate(
      flag_no_spots = is.na(n_spots) | n_spots == 0,
      
      flag_lt100 = !flag_no_spots & n_spots < config$sensitivity_min_spots,
      
      flag_100_199 = !flag_no_spots & n_spots >= config$sensitivity_min_spots &
        n_spots < config$primary_min_spots,
      
      eligible_primary = !flag_no_spots & n_spots >= config$primary_min_spots &
        pct_classification_complete == 100,
      
      # All candidate predictors use the same spatial eligibility.
      # Module coverage is reported separately and is not used as a subarray exclusion rule in the primary validation.
      eligible_tumor = eligible_primary,
      eligible_lymphocyte = eligible_primary,
      eligible_tls = eligible_primary,
      
      review_status = dplyr::case_when(
        flag_no_spots ~ "No retained spots",
        
        subarray_id == "CN26_D1" ~ "Exclude: only 8 retained spots",
        
        flag_lt100 ~ "Exclude from primary spatial analysis: < 100 spots",
        
        flag_100_199 ~ "Sensitivity only: 100-199 spots",
        
        !eligible_primary ~ "Review: incomplete classification or other QC issue",
        
        TRUE ~ "Primary eligible"
      )
    )
}

build_validation_data <- function(annotated_spots, eligibility) {
  required <- c(
    "subarray_id", "count_id", 
    "manual_tumor", "manual_lymphocyte", "manual_tls", 
    "manual_tumor_positive","manual_lymphocyte_positive", "manual_tls_positive",
    "hist_score_Tumor", "hist_score_Lymphocyte", "hist_score_Lymphoid.nodule",
    "score_TUMOR_EPITHELIAL", "score_TUMOR_EPITHELIAL_z", "score_TNBC_BASAL",
    "score_TNBC_BASAL_z", "score_T_CELL", "score_T_CELL_z","score_B_CELL",
    "score_B_CELL_z", "score_TLS_WANG_30", "score_TLS_WANG_30_z","score_TLS_LIKE",
    "score_TLS_LIKE_z", "score_PLASMA", "score_PLASMA_z"
  )

  missing <- setdiff(required, names(annotated_spots))
  if (length(missing) > 0) {
    stop(
      "annotated_spots is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  annotated_spots |>
    dplyr::left_join(
      eligibility |>
        dplyr::select(subarray_id, eligible_primary, eligible_tumor,
          eligible_lymphocyte, eligible_tls, review_status),
      by = "subarray_id") |>
    dplyr::group_by(subarray_id) |>
    dplyr::mutate(
      sample_id = as.character(count_id),
      
      tumor_expr_max_z = pmax(
        score_TUMOR_EPITHELIAL_z, score_TNBC_BASAL_z, na.rm = TRUE),
      
      lymphocyte_expr_broad_z = rowMeans(
        cbind(score_T_CELL_z, score_B_CELL_z), na.rm = TRUE),
      
      # Ranks must be calculated within each subarray.
      tumor_expr_rank = pmax(
        rank01(score_TUMOR_EPITHELIAL), rank01(score_TNBC_BASAL), na.rm = TRUE),
      
      lymphocyte_expr_broad_rank = rowMeans(
        cbind(rank01(score_T_CELL), rank01(score_B_CELL)), na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      tumor_expr_max_z =
        ifelse(is.finite(tumor_expr_max_z), tumor_expr_max_z, NA_real_),
      immune_expr_broad_z =
        ifelse(is.finite(lymphocyte_expr_broad_z), lymphocyte_expr_broad_z, NA_real_),
      tumor_expr_rank =
        ifelse(is.finite(tumor_expr_rank), tumor_expr_rank, NA_real_),
      immune_expr_broad_rank =
        ifelse(is.finite(lymphocyte_expr_broad_rank), lymphocyte_expr_broad_rank, NA_real_)
    )
}

PREDICTOR_SETS <- list(
  tumor = c(
    "hist_score_Tumor",
    "score_TUMOR_EPITHELIAL",
    "score_TUMOR_EPITHELIAL_z",
    "score_TNBC_BASAL",
    "score_TNBC_BASAL_z",
    "tumor_expr_max_z",
    "tumor_expr_rank"
  ),

  lymphocyte = c(
    "hist_score_Lymphocyte",
    "score_T_CELL",
    "score_T_CELL_z",
    "score_B_CELL",
    "score_B_CELL_z",
    "lymphocyte_expr_broad_z",
    "lymphocyte_expr_broad_rank"
  ),

  tls = c(
    "hist_score_Lymphoid.nodule",
    "score_TLS_WANG_30",
    "score_TLS_WANG_30_z",
    "score_TLS_LIKE",
    "score_TLS_LIKE_z",
    "score_B_CELL",
    "score_B_CELL_z",
    "score_PLASMA",
    "score_PLASMA_z"
  )
)

OUTCOME_MAP <- list(
  tumor = list(
    fraction = "manual_tumor",
    binary = "manual_tumor_positive",
    eligibility = "eligible_tumor"
  ),

  lymphocyte = list(
    fraction = "manual_lymphocyte",
    binary = "manual_lymphocyte_positive",
    eligibility = "eligible_lymphocyte"
  ),

  tls = list(
    fraction = "manual_tls",
    binary = "manual_tls_positive",
    eligibility = "eligible_tls"
  )
)

evaluate_one_predictor <- function(dat, target, predictor, config = CONFIG) {
  info <- OUTCOME_MAP[[target]]
  
  dat |>
    dplyr::filter(.data[[info$eligibility]]) |>
    dplyr::group_by(sample_id, subarray_id) |>
    dplyr::group_modify(function(df, key) {
      x <- df[[predictor]]
      y <- df[[info$fraction]]
      y_bin <- df[[info$binary]]
      
      data.frame(
        target = target,
        predictor = predictor,
        n_spots = nrow(df),
        n_positive = sum(y_bin, na.rm = TRUE),
        n_negative = sum(!y_bin, na.rm = TRUE),
        rho = safe_cor(x,y, config$min_complete_spots_for_rho),
        auc = safe_auc(x, y_bin, config$min_positive_for_auc, config$min_negative_for_auc)
      )
    }) |>
    dplyr::ungroup()
}

evaluate_all_predictors <- function(dat, config = CONFIG) {
  out <- list()
  
  for (target in names(PREDICTOR_SETS)) {
    for (predictor in PREDICTOR_SETS[[target]]) {
      if (!predictor %in% names(dat)) next
      
      out[[paste(target, predictor, sep = "__")]] <-
        evaluate_one_predictor(dat, target, predictor, config)
    }
  }
  
  data.table::rbindlist(out, fill = TRUE) |>
    as.data.frame()
}

summarize_performance <- function(perf) {
  perf |>
    dplyr::group_by(target, predictor) |>
    dplyr::summarise(
      n_subarrays = dplyr::n_distinct(subarray_id),
      n_evaluable_rho = sum(is.finite(rho)),
      n_evaluable_auc = sum(is.finite(auc)),
      median_rho = median(rho, na.rm = TRUE),
      q25_rho = safe_quantile(rho, 0.25),
      q75_rho = safe_quantile(rho, 0.75),
      pct_positive_rho = mean(rho > 0, na.rm = TRUE) * 100,
      median_auc = median(auc, na.rm = TRUE),
      q25_auc = safe_quantile(auc, 0.25),
      q75_auc = safe_quantile(auc, 0.75),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      target,
      dplyr::desc(median_auc),
      dplyr::desc(median_rho)
    )
}

make_summary_plot <- function(summary_tbl, metric, y_label) {
  ggplot2::ggplot(
    summary_tbl,
    ggplot2::aes(
      x = reorder(predictor, .data[[metric]]),
      y = .data[[metric]]
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~target, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = NULL,
      y = y_label,
      title = paste("Candidate marker-score validation:", y_label)
    )
}

run_validation <- function(config = CONFIG) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  setwd(config$project_dir)
  output_dir <- config$working_output_dir
  results_dir <- config$validation_results_dir

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  qc <- data.table::fread(
    file.path(output_dir, "cohort_qc_summary.csv")
  ) |> as.data.frame()

  coverage <- data.table::fread(
    file.path(output_dir, "cohort_module_coverage.csv")
  ) |> as.data.frame()

  annotated_spots <- build_annotated_spots(config)
  eligibility <- build_eligibility_table(qc, coverage, config)

  data.table::fwrite(
    eligibility,
    file.path(results_dir, "subarray_eligibility_for_review.csv")
  )

  validation_data <- build_validation_data(annotated_spots, eligibility)
  saveRDS(
    validation_data,
    file.path(results_dir, "annotated_validation_spots.rds"),
    compress = TRUE
  )

  performance_by_subarray <- evaluate_all_predictors(validation_data, config)
  data.table::fwrite(
    performance_by_subarray,
    file.path(results_dir, "predictor_performance_by_subarray.csv")
  )

  performance_summary <- summarize_performance(performance_by_subarray)
  data.table::fwrite(
    performance_summary,
    file.path(results_dir, "predictor_performance_summary.csv")
  )

  provisional_best <- performance_summary |>
    dplyr::group_by(target) |>
    dplyr::arrange(
      dplyr::desc(median_auc),
      dplyr::desc(median_rho),
      dplyr::desc(n_evaluable_auc)
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  data.table::fwrite(
    provisional_best,
    file.path(results_dir, "provisional_best_predictors.csv")
  )

  plot_auc <- make_summary_plot(
    performance_summary, "median_auc", "Median subarray-level AUROC"
  )

  plot_rho <- make_summary_plot(
    performance_summary, "median_rho", "Median subarray-level Spearman correlation"
  )

  ggplot2::ggsave(
    file.path(results_dir, "candidate_predictor_median_auc.pdf"),
    plot_auc, width = 10, height = 7
  )

  ggplot2::ggsave(
    file.path(results_dir, "candidate_predictor_median_rho.pdf"),
    plot_rho, width = 10, height = 7
  )

  cohort_summary <- data.frame(
    metric = c(
      "Processed subarrays",
      "Primary eligible subarrays",
      "Sensitivity-only subarrays",
      "Subarrays with <100 spots",
      "Annotated subarrays loaded",
      "Annotated primary-eligible subarrays",
      "Annotated source samples loaded",
      "Annotated primary-eligible source samples"
    ),

    value = c(nrow(qc),
              sum(eligibility$eligible_primary, na.rm = TRUE),
              sum(eligibility$flag_100_199, na.rm = TRUE),
              sum(eligibility$flag_lt100, na.rm = TRUE),
              dplyr::n_distinct(validation_data$subarray_id),
              dplyr::n_distinct(validation_data$subarray_id[validation_data$eligible_primary]),
              dplyr::n_distinct(validation_data$sample_id),
              dplyr::n_distinct(validation_data$sample_id[validation_data$eligible_primary])
    )
  )

  data.table::fwrite(
    cohort_summary,
    file.path(results_dir, "cohort_validation_summary.csv")
  )

  invisible(
    list(
      eligibility = eligibility,
      validation_data = validation_data,
      performance_by_subarray = performance_by_subarray,
      performance_summary = performance_summary,
      provisional_best = provisional_best,
      cohort_summary = cohort_summary
    )
  )
}

if (sys.nframe() == 0) {
  results <- run_validation(CONFIG)
}
