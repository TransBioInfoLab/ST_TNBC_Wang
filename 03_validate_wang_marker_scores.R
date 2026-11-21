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
  x <- as.data.frame(readRDS(file), check.names = FALSE)
  
  sid <- sub("_spots\\.rds$", "", basename(file))
  
  if (!"spot_id" %in% names(x) ||
      anyNA(x$spot_id) ||
      !all(startsWith(
        as.character(x$spot_id), paste0(sid, ".")
      ))) {
    stop("Spot IDs do not match the array filename: ", file)
  }
  
  reconstructed <- !"subarray_id" %in% names(x)
  if (reconstructed) x$subarray_id <- sid
  
  if (!all(x$subarray_id %in% sid)) {
    stop("Wrong subarray_id in ", file)
  }
  
  x$subarray_id_reconstructed <- reconstructed
  
  fixed <- c(
    "spot_id", "run_id", "count_id", "subarray_id",
    "subarray_id_reconstructed", "subject_id",
    "specimen_type", "baseline_primary_specimen",
    "qc_pass", "has_manual_annotation",
    "module_score_version"
  )
  
  required <- c(
    fixed,
    "manual_annotated_pixels",
    "manual_annotation_supported",
    "manual_annotation_min_pixels"
  )
  
  if (!all(required %in% names(x))) {
    stop(
      "Rebuild this spot table with corrected Step 02: ",
      file
    )
  }
  
  keep <- unique(c(
    fixed,
    grep("^(hist_|manual_|score_)", names(x), value = TRUE)
  ))
  
  out <- x[, keep, drop = FALSE]
  
  scalar <- vapply(out, function(z) {
    is.atomic(z) &&
      is.null(dim(z)) &&
      length(z) == nrow(out)
  }, logical(1))
  
  if (any(!scalar) || anyDuplicated(out$spot_id)) {
    stop("Invalid spot table: ", file)
  }
  
  out
}

build_annotated_spots <- function(qc, config) {
  ids <- as.character(qc$subarray_id[qc$has_manual_annotation %in% TRUE])
  
  files <- file.path(config$working_output_dir, "spot_tables", paste0(ids, "_spots.rds"))
  
  if (!length(files) || any(!file.exists(files))) {
    stop("Annotated spot tables are absent or incomplete.")
  }
  
  x <- as.data.frame(data.table::rbindlist(
    lapply(files, read_validation_spot_table),
    fill = TRUE
  ))
  
  if (anyDuplicated(x$spot_id) || !all(x$run_id %in% config$run_id)) {
    stop(
      "Annotated spots contain duplicate IDs ",
      "or a different run_id."
    )
  }
  
  observed <- table(x$subarray_id)
  expected <- qc$n_spots[match(names(observed), qc$subarray_id)]
  
  if (!setequal(names(observed), ids) ||
      any(as.integer(observed) != expected)) {
    stop("Annotated spot counts do not match the processing QC.")
  }
  
  supported <- is.finite(x$manual_annotated_pixels) &
    x$manual_annotated_pixels > config$min_annotation_biological_pixels
  
  if (anyNA(x$manual_annotation_supported) ||
      !all(x$manual_annotation_supported == supported) ||
      !all(x$manual_annotation_min_pixels %in%
           config$min_annotation_biological_pixels)) {
    stop(
      "Annotation support does not match ",
      "the frozen processing config."
    )
  }
  
  x$manual_tumor <- ifelse( 
    supported, 
    x$manual_fraction_Tumor + x$manual_fraction_Tumor.region,
    NA_real_
  )
  
  x$manual_lymphocyte <- ifelse(
    supported,
    x$manual_fraction_Lymphocyte + 0.5 * x$manual_fraction_High.TIL.stroma,
    NA_real_
  )
  
  x$manual_tls <- ifelse(
    supported,
    x$manual_fraction_Lymphoid.nodule,
    NA_real_
  )
  
  for (target in c("tumor", "lymphocyte", "tls")) {
    x[[paste0("manual_", target, "_positive")]] <-
      x[[paste0("manual_", target)]] > config$manual_positive_threshold
  }
  x
}

attach_source_loo <- function(dat, config) {
  ids <- readRDS(config$ids_file)
  k <- which(ids$hasAnnot %in% TRUE)
  
  count_id <- as.character(ids$id[k])
  subarray_id <- rownames(ids)[k]
  
  if (anyDuplicated(count_id)) {
    stop(
      "Source LOO reconstruction assumes ",
      "one annotated array per specimen."
    )
  }
  
  files <- file.path(config$annotation_dir, paste0("TNBC", count_id, ".RDS"))
  annotations <- lapply(files, function(f) readRDS(f)$annots)
  n <- do.call(rbind, annotations)
  
  spot_id <- unlist(Map(function(a, sid) {
    paste0(sid, ".", rownames(a))
  }, annotations, subarray_id), use.names = FALSE)
  
  specimens <- rep(count_id, vapply(annotations, nrow, integer(1)))
  n <- cbind(n, Stroma = NA_real_)
  
  n <- n[, !colnames(n) %in% c("Nothing", "Artefacts", "Hole (whitespace)"), drop = FALSE]
  n[, "Tumor"] <- n[, "Tumor"] + n[, "Tumor region"]
  n[, "Stroma"] <- n[, "Stroma cell"] + n[, "Low TIL stroma"] +
    n[, "High TIL stroma"] / 2 + n[, "Acellular stroma"]
  n[, "Lymphocyte"] <- n[, "Lymphocyte"] + n[, "High TIL stroma"] / 2
  
  n <- n[, !colnames(n) %in% c("Low TIL stroma", "High TIL stroma", "Stroma cell", 
                               "Acellular stroma", "Tumor region"), drop = FALSE]
  
  # Fixed SOURCE-model filter.
  valid <- rowSums(n) > 1000
  n <- n[valid, , drop = FALSE]
  n <- n / rowSums(n)
  spot_id <- spot_id[valid]
  specimens <- specimens[valid]
  
  loo <- new.env()
  load(config$heldout_file, envir = loo)
  
  if (!identical(dim(n), dim(loo$n)) ||
      !identical(dimnames(n), dimnames(loo$n)) ||
      any(!is.finite(n)) ||
      any(!is.finite(loo$n)) ||
      max(abs(n - loo$n)) > 1e-12 ||
      anyDuplicated(spot_id)) {
    stop("Cannot establish the source LOO row ordering.")
  }
  
  m <- match(dat$spot_id, spot_id)
  source_supported <- is.finite(dat$manual_annotated_pixels) & dat$manual_annotated_pixels > 1000
  
  if (any(source_supported & is.na(m)) ||
      any(
        as.character(dat$count_id[!is.na(m)]) !=
        specimens[m[!is.na(m)]]
      )) {
    stop(
      "Source LOO identifiers do not match ",
      "the processed spots."
    )
  }
  
  targets <- c(tumor = "Tumor", lymphocyte = "Lymphocyte", tls = "Lymphoid nodule")
  
  for (target in names(targets)) {
    class <- targets[[target]]
    prediction <- loo$loos[[class]]
    
    if (!is.numeric(prediction) ||length(prediction) != nrow(n)) {
      stop("Invalid source LOO vector: ", class)
    }
    
    dat[[paste0("hist_loo_", make.names(class))]] <- prediction[m]
    truth <- n[m, class]
    observed <- dat[[paste0("manual_", target)]]
    ok <- is.finite(truth) & is.finite(observed)
    
    if (any(abs(truth[ok] - observed[ok]) > 1e-12)) {
      stop(
        "Manual target differs from the source definition: ",
        target
      )
    }
  }
  
  list(
    data = dat,
    source_files = c(config$ids_file, files, config$heldout_file)
  )
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
  
  if (sum(y == TRUE) < min_pos || sum(y == FALSE) < min_neg) {
    return(NA_real_)
  }
  
  if (length(unique(x)) == 1L) return(0.5)
  
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(
      y, x, direction = "<", quiet = TRUE
    ))),
    error = function(e) NA_real_
  )
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

build_eligibility_table <- function(qc, coverage, config, status) {
  ledger <- as.data.frame(status[, c(
    "run_id", "count_id", "subarray_id", "success"
  ), with = FALSE])
  
  names(ledger)[names(ledger) == "success"] <- "processing_success"
  q <- qc[, setdiff(names(qc), c("run_id", "count_id")), drop = FALSE]
  x <- dplyr::left_join(ledger, q, by = "subarray_id")
  
  # Include hierarchy for the known empty array.
  missing <- is.na(x$subject_id)
  x$subject_id[missing] <- as.character(x$count_id[missing])
  x$subject_id[x$count_id == "30"] <- "53"
  x$subject_id[x$count_id == "66"] <- "58"
  
  x$specimen_type <- ifelse(
    x$count_id %in% c("30", "66"), "relapse", "primary"
  )
  x$baseline_primary_specimen <- x$specimen_type == "primary"
  
  coverage_wide <- coverage |>
    dplyr::select(
      subarray_id, module,
      fraction_available, fraction_detected,
      score_eligible
    ) |>
    tidyr::pivot_wider(
      names_from = module,
      values_from = c(
        fraction_available,
        fraction_detected,
        score_eligible
      )
    )
  
  dplyr::left_join(
    x, coverage_wide, by = "subarray_id"
  ) |>
    dplyr::mutate(
      flag_no_spots = is.na(n_spots) | n_spots == 0,
      
      flag_lt100 = !flag_no_spots &
        n_spots < config$sensitivity_min_spots,
      
      flag_100_199 = !flag_no_spots &
        n_spots >= config$sensitivity_min_spots &
        n_spots < config$primary_min_spots,
      
      eligible_spatial_main =
        processing_success %in% TRUE &
        !flag_no_spots &
        n_spots >= config$primary_min_spots &
        pct_classification_complete %in% 100,
      
      # Compatibility alias: analysis tier, not specimen type.
      eligible_primary = eligible_spatial_main,
      
      eligible_baseline_primary = eligible_spatial_main & baseline_primary_specimen,
      
      eligible_tumor = eligible_baseline_primary,
      eligible_lymphocyte = eligible_baseline_primary,
      eligible_tls = eligible_baseline_primary,
      
      review_status = dplyr::case_when(
        flag_no_spots ~ "No retained spots",
        flag_lt100 ~ "Below sensitivity spot threshold",
        flag_100_199 ~ "Sensitivity only",
        !eligible_spatial_main ~
          "Not spatial-main eligible",
        !baseline_primary_specimen ~
          "Relapse; excluded from primary validation",
        TRUE ~ "Primary specimen eligible"
      )
    )
}

build_validation_data <- function(
    annotated_spots, eligibility
) {
  required <- c(
    "subject_id", "specimen_type", "qc_pass",
    "manual_annotation_supported",
    "hist_score_Tumor",
    "hist_score_Lymphocyte",
    "hist_score_Lymphoid.nodule",
    "hist_loo_Tumor",
    "hist_loo_Lymphocyte",
    "hist_loo_Lymphoid.nodule",
    "score_TUMOR_EPITHELIAL",
    "score_TUMOR_EPITHELIAL_z",
    "score_TNBC_BASAL", "score_TNBC_BASAL_z",
    "score_T_CELL", "score_T_CELL_z",
    "score_B_CELL", "score_B_CELL_z",
    "score_TLS_WANG_30", "score_TLS_WANG_30_z",
    "score_TLS_LIKE", "score_TLS_LIKE_z",
    "score_PLASMA", "score_PLASMA_z"
  )
  
  if (!all(required %in% names(annotated_spots))) {
    stop("Validation predictors/metadata are missing.")
  }
  
  annotated_spots |>
    dplyr::left_join(
      eligibility |>
        dplyr::select(
          subarray_id,
          eligible_primary,
          eligible_spatial_main,
          eligible_baseline_primary,
          eligible_tumor,
          eligible_lymphocyte,
          eligible_tls,
          review_status
        ),
      by = "subarray_id"
    ) |>
    dplyr::group_by(subarray_id) |>
    dplyr::mutate(
      sample_id = as.character(count_id),
      
      # Require all components of each composite.
      tumor_expr_max_z = pmax(
        score_TUMOR_EPITHELIAL_z,
        score_TNBC_BASAL_z
      ),
      
      lymphocyte_expr_broad_z = rowMeans(cbind(
        score_T_CELL_z, score_B_CELL_z
      )),
      
      tumor_expr_rank = pmax(
        rank01(score_TUMOR_EPITHELIAL),
        rank01(score_TNBC_BASAL)
      ),
      
      lymphocyte_expr_broad_rank = rowMeans(cbind(
        rank01(score_T_CELL),
        rank01(score_B_CELL)
      ))
    ) |>
    dplyr::ungroup()
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

PREDICTOR_SETS$tumor <- c("hist_loo_Tumor", PREDICTOR_SETS$tumor)

PREDICTOR_SETS$lymphocyte <- c("hist_loo_Lymphocyte", PREDICTOR_SETS$lymphocyte)

PREDICTOR_SETS$tls <- c("hist_loo_Lymphoid.nodule", PREDICTOR_SETS$tls)

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
  predictors <- PREDICTOR_SETS[[target]]
  
  type <- if (startsWith(predictor, "hist_loo_")) {
    "source_specimen_LOO"
  } else if (startsWith(predictor, "hist_score_")) {
    "source_full_fit"
  } else {
    "fixed_module_annotation_agreement"
  }
  
  dat |>
    dplyr::filter(
      .data[[info$eligibility]] %in% TRUE
    ) |>
    dplyr::group_by(
      subject_id, sample_id, subarray_id
    ) |>
    dplyr::group_modify(function(df, key) {
      x <- df[[predictor]]
      y <- df[[info$fraction]]
      y_bin <- df[[info$binary]]
      
      eligible <-
        df$qc_pass %in% TRUE &
        df$manual_annotation_supported %in% TRUE &
        is.finite(y) & !is.na(y_bin)
      
      common <- eligible & Reduce(
        "&",
        lapply(
          df[, predictors, drop = FALSE], is.finite)
      )
      
      data.table::rbindlist(lapply(
        c("common_predictors", "available_predictor"),
        function(set) {
          keep <- if (set == "common_predictors") {
            common
          } else {
            eligible & is.finite(x)
          }
          
          data.frame(
            target = target,
            predictor = predictor,
            prediction_type = type,
            comparison_set = set,
            
            linked_specimen_in_training = type == "source_specimen_LOO" &
              key$subject_id %in% c("53", "58"),
            
            n_target_eligible = sum(eligible),
            n_spots = sum(keep),
            n_positive = sum(y_bin[keep]),
            n_negative = sum(!y_bin[keep]),
            
            rho = safe_cor(x[keep], y[keep], config$min_complete_spots_for_rho),
            
            auc = safe_auc(x[keep], y_bin[keep],
              config$min_positive_for_auc, config$min_negative_for_auc)
          )
        }
      )) |>
        as.data.frame()
    }) |>
    dplyr::ungroup()
}

evaluate_all_predictors <- function(dat, config = CONFIG) {
  out <- list()
  
  for (target in names(PREDICTOR_SETS)) {
    for (predictor in PREDICTOR_SETS[[target]]) {
      if (!predictor %in% names(dat)) {
        stop("Missing predictor: ", predictor)
      }
      
      out[[paste(target, predictor, sep = "__")]] <-
        evaluate_one_predictor(dat, target, predictor, config)
    }
  }
  
  data.table::rbindlist(out, fill = TRUE) |>
    as.data.frame()
}

summarize_performance <- function(perf) {
  perf |>
    dplyr::group_by(
      target, comparison_set, subarray_id
    ) |>
    dplyr::mutate(
      matched_rho = all(is.finite(rho)),
      matched_auc = all(is.finite(auc))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      rho_used = ifelse(
        comparison_set == "common_predictors" & !matched_rho,
        NA_real_, rho
      ),
      auc_used = ifelse(
        comparison_set == "common_predictors" & !matched_auc,
        NA_real_, auc
      )
    ) |>
    dplyr::group_by(target, comparison_set, prediction_type, predictor) |>
    dplyr::summarise(
      n_subarrays = dplyr::n_distinct(subarray_id),
      n_patients = dplyr::n_distinct(subject_id),
      
      n_evaluable_rho = sum(is.finite(rho_used)),
      n_evaluable_auc = sum(is.finite(auc_used)),
      
      median_rho = safe_quantile(rho_used, 0.5),
      q25_rho = safe_quantile(rho_used, 0.25),
      q75_rho = safe_quantile(rho_used, 0.75),
      
      median_auc = safe_quantile(auc_used, 0.5),
      q25_auc = safe_quantile(auc_used, 0.25),
      q75_auc = safe_quantile(auc_used, 0.75),
      
      .groups = "drop"
    )
}

make_summary_plot <- function(summary_tbl, metric, y_label) {
  summary_tbl <- dplyr::filter(summary_tbl, comparison_set == "common_predictors")
  
  ggplot2::ggplot(
    summary_tbl,
    ggplot2::aes(
      x = reorder(predictor, .data[[metric]]),
      y = .data[[metric]],
      fill = prediction_type
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      ~target, scales = "free_y"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      x = NULL,
      y = y_label,
      title = paste(
        "Primary specimens;",
        "matched spots and evaluable sections"
      ),
      fill = "Prediction source"
    )
}

run_validation <- function(config = CONFIG) {
  used_config_file <- normalizePath(config_file, mustWork = TRUE)
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(config$project_dir)
  
  root <- normalizePath(config$working_output_dir, mustWork = TRUE)
  
  results_dir <- file.path(
    normalizePath(
      dirname(config$validation_results_dir),
      mustWork = TRUE
    ),
    basename(config$validation_results_dir)
  )
  
  if (!startsWith(results_dir, paste0(root, "/")) ||
      (dir.exists(results_dir) &&
       length(list.files(results_dir, all.files = TRUE, no.. = TRUE)))) {
    stop(
      "Choose an empty validation directory ",
      "inside this processing release."
    )
  }
  
  logs <- file.path(root, "logs")
  proc <- readRDS(file.path(logs, "processing_config.rds"))
  
  report <- data.table::fread(file.path(logs, "processing_run_summary.csv"))
  
  if (!all(c("run_id", "run_status") %in% names(report)) ||
  nrow(report) != 1L || !all(report$run_status %in% "complete") || 
  !all(report$run_id %in% proc$run_id)) {
    stop("Step 02 must have completed successfully.")
  }
  
  if (!identical(as.numeric(config$min_annotation_biological_pixels),
    as.numeric(proc$min_annotation_biological_pixels))) {
    stop(
      "Use the same annotation-support threshold ",
      "as Step 02."
    )
  }
  
  config$run_id <- proc$run_id
  config$working_output_dir <- root
  config$ids_file <- proc$ids_file
  config$annotation_dir <- proc$annotation_dir
  config$validation_results_dir <- results_dir
  config$annotated_spots_file <- file.path(results_dir, "annotated_spots.rds")
  
  config$validation_run_id <- paste0(
    "validation_", config$run_id, "_",
    format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
    "Z"
  )
  
  qc_file <- file.path(root, "cohort_qc_summary.csv")
  coverage_file <- file.path(root, "cohort_module_coverage.csv")
  status_file <- file.path(
    root, "status",
    "subarray_processing_status_latest.csv"
  )
  
  qc <- as.data.frame(data.table::fread(qc_file))
  coverage <- as.data.frame(data.table::fread(coverage_file))
  status <- data.table::fread(status_file)
  
  for (x in list(qc, coverage, status)) {
    if (!"run_id" %in% names(x) ||
        !all(x$run_id %in% config$run_id)) {
      stop(
        "Processing tables contain missing ",
        "or inconsistent run identifiers."
      )
    }
  }
  
  good <- status[success %in% TRUE, subarray_id]
  defs <- readRDS(file.path(logs, "scoring_definitions.rds"))
  
  if (anyDuplicated(qc$subarray_id) ||
      !setequal(good, qc$subarray_id) ||
      nrow(coverage) !=
      length(good) * length(defs$gene_sets) ||
      !all(coverage$subarray_id %in% good) ||
      !all(coverage$module %in% names(defs$gene_sets)) ||
      anyDuplicated(
        coverage[, c("subarray_id", "module")]
      )) {
    stop(
      "QC/coverage records do not match ",
      "the completed processing run."
    )
  }
  
  primary_annotated <- qc[qc$has_manual_annotation %in% TRUE & qc$specimen_type == "primary",]
  
  if (anyDuplicated(primary_annotated$subject_id)) {
    stop(
      "This validation expects one annotated ",
      "primary section per patient."
    )
  }
  
  annotated_spots <- build_annotated_spots(qc, config)
  loo <- attach_source_loo(annotated_spots, config)
  eligibility <- build_eligibility_table(qc, coverage, config, status)
  validation_data <- build_validation_data(loo$data, eligibility)
  
  # Hash inputs actually consumed by validation.
  script_file <- normalizePath(
    "03_validate_wang_marker_scores.R",
    mustWork = TRUE
  )
  
  spot_files <- file.path(
    root, "spot_tables",
    paste0(
      qc$subarray_id[qc$has_manual_annotation %in% TRUE],
      "_spots.rds"
    )
  )
  
  paths <- unique(normalizePath(
    c(
      used_config_file,
      script_file,
      qc_file,
      coverage_file,
      status_file,
      file.path(logs, c(
        "processing_config.rds",
        "processing_run_summary.csv",
        "scoring_definitions.rds",
        "processing_input_manifest.csv"
      )),
      loo$source_files,
      spot_files
    ),
    mustWork = TRUE
  ))
  
  inputs <- data.frame(path = paths, 
                       md5 = unname(tools::md5sum(paths)),
                       stringsAsFactors = FALSE
  )
  
  if (anyNA(inputs$md5)) {
    stop("Could not checksum validation inputs.")
  }
  
  previous <- data.table::fread(file.path(logs, "processing_input_manifest.csv"))
  previous <- previous[
    !input %in% c("config", "runner", "helpers", "gene_sets")
  ]
  
  j <- match(inputs$path, previous$path)
  if (any(!is.na(j) & inputs$md5 != previous$md5[j])) {
    stop(
      "A consumed source-data file ",
      "changed since Step 02."
    )
  }
  
  performance <- evaluate_all_predictors(validation_data, config)
  summary <- summarize_performance(performance)
  
  flow <- data.frame(
    metric = c(
      "Processed subarrays",
      "Independent processed patients",
      "Spatial-main eligible subarrays",
      "Baseline-primary eligible subarrays",
      "Annotated specimens",
      "Annotated primary patients",
      "Supported annotation spots",
      "Spots with source tumor LOO predictions"
    ),
    value = c(
      nrow(qc),
      dplyr::n_distinct(qc$subject_id),
      sum(eligibility$eligible_spatial_main),
      sum(eligibility$eligible_baseline_primary),
      dplyr::n_distinct(annotated_spots$count_id),
      dplyr::n_distinct(
        annotated_spots$subject_id[
          annotated_spots$specimen_type == "primary"
        ]
      ),
      sum(annotated_spots$manual_annotation_supported),
      sum(is.finite(loo$data$hist_loo_Tumor))
    )
  )
  
  stamp <- function(x) {
    x <- as.data.frame(x)
    x$run_id <- config$run_id
    x$validation_run_id <- config$validation_run_id
    x
  }
  
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  tables <- list(
    subarray_eligibility_for_review = eligibility,
    predictor_performance_by_subarray = performance,
    predictor_performance_summary = summary,
    cohort_validation_summary = flow,
    validation_input_manifest = inputs
  )
  
  for (nm in names(tables)) {
    data.table::fwrite(
      stamp(tables[[nm]]),
      file.path(results_dir, paste0(nm, ".csv"))
    )
  }
  
  saveRDS(stamp(annotated_spots), config$annotated_spots_file, compress = TRUE)
  
  saveRDS(
    stamp(validation_data),
    file.path(
      results_dir,
      "annotated_validation_spots.rds"
    ),
    compress = TRUE
  )
  
  saveRDS(
    config,
    file.path(results_dir, "validation_config.rds")
  )
  
  if (!file.copy(
    script_file,
    file.path(results_dir, "validation_script.R"),
    overwrite = FALSE
  )) {
    stop("Could not save the validation script snapshot.")
  }
  
  capture.output(
    sessionInfo(),
    file = file.path(results_dir, "validation_session_info.txt")
  )
  
  for (metric in c("auc", "rho")) {
    plot <- make_summary_plot(
      summary, paste0("median_", metric),
      if (metric == "auc") {
        "Median primary-patient AUROC"
      } else {
        "Median primary-patient Spearman correlation"
      }
    )
    
    ggplot2::ggsave(
      file.path(results_dir, paste0("candidate_predictor_median_", metric, ".pdf")),
      plot, width = 10, height = 7
    )
  }
  
  # Completion record is written last.
  data.table::fwrite(
    stamp(data.frame(
      validation_status = "complete",
      main_comparison = "primary specimens; matched supported spots",
      heldout_unit = "source specimen, not linked patient"
    )),
    file.path(results_dir, "validation_run_summary.csv")
  )
  
  invisible(list(
    eligibility = eligibility,
    validation_data = validation_data,
    performance_by_subarray = performance,
    performance_summary = summary,
    cohort_summary = flow
  ))
}

if (sys.nframe() == 0) {
  results <- run_validation(CONFIG)
}

