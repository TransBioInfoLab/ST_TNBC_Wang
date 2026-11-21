#!/usr/bin/env Rscript
# Step 05: final release assembly after Step 04 overlay QA and recorded review.
# Usage: Rscript 05_prepare_wang_release.R [--dry-run] [--overlay-dir=QA_RUN]
#        [--processing-dir=PATH] [--validation-dir=PATH] [--release-dir=PATH]
# Default final folder: outputs_wang_v1.1/Wang_TNBC_v1.1.0 (same existing layout).
# Processing/validation inputs are read only. The release is built in staging
# and renamed only after verification. An existing final folder is preserved
# in a sibling backup before the staged folder takes its original name.
# Pending alignment/orientation decisions block final assembly; dry-run writes nothing.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

require_release <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

find_wang_config <- function() {
  configured <- Sys.getenv("WANG_TNBC_CONFIG", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  folder <- if (length(args)) {
    dirname(normalizePath(sub("^--file=", "", args[1]), mustWork = TRUE))
  } else getwd()
  normalizePath(file.path(folder, "00_wang_config.R"), mustWork = TRUE)
}

canonical_release_version <- function(x) {
  x <- sub("^v", "", as.character(x))
  require_release(length(x) == 1L && !is.na(x), "Supply one release version.")
  if (grepl("^[0-9]+\\.[0-9]+$", x)) x <- paste0(x, ".0")
  require_release(grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", x),
                  "Use a release version such as 1.1.0.")
  x
}

relative_release_path <- function(x) {
  require_release(length(x) == 1L && !is.na(x) && nzchar(x),
                  "Missing release file path.")
  x <- gsub("\\\\", "/", x)
  require_release(!grepl("^/|^[A-Za-z]:|(^|/)\\.\\.?(/|$)", x),
                  paste("Unsafe release path:", x))
  x
}

owned_file <- function(root, relative) {
  relative <- relative_release_path(relative)
  path <- normalizePath(file.path(root, relative), winslash = "/", mustWork = TRUE)
  require_release(startsWith(path, paste0(root, "/")),
                  paste("File escapes processing directory:", relative))
  path
}

file_md5 <- function(paths) {
  hashes <- unname(tools::md5sum(paths))
  require_release(!anyNA(hashes), "Could not checksum every required file.")
  hashes
}

project_path <- function(path, project) {
  if (!grepl("^/", path)) path <- file.path(project, path)
  path
}

resolve_processing_root <- function(config, version) {
  if (!is.null(config$processing_input_dir)) {
    return(normalizePath(project_path(config$processing_input_dir, config$project_dir), mustWork = TRUE))
  }
  working <- project_path(config$working_output_dir, config$project_dir)
  candidates <- c(working, file.path(working, paste0("Wang_TNBC_v", version)))
  found <- candidates[file.exists(file.path(candidates, "logs", "processing_config.rds"))]
  require_release(length(found) > 0L, "Completed Step 02 output not found; use --processing-dir=PATH.")
  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

read_frozen_definitions <- function(root, good, vinputs) {
  relative <- "logs/scoring_definitions.rds"
  recorded <- vinputs[basename(path) == "scoring_definitions.rds"]
  require_release(nrow(recorded) == 1L, "Step 03 did not record frozen scoring definitions.")
  source <- file.path(root, relative)
  recovered <- !file.exists(source)
  evidence <- relative
  if (recovered) {
    # The previous assembler omitted this file. Recover only if its original
    # BYTE checksum is reproduced exactly; do not create new scoring definitions.
    evidence <- file.path("module_coverage", paste0(good$subarray_id[1], "_scoring_details.rds"))
    info <- readRDS(owned_file(root, evidence))$module_scoring
    source <- tempfile("wang_frozen_definitions_", fileext = ".rds")
    saveRDS(list(gene_sets = info$gene_sets, alias_groups = info$alias_groups), source)
  }
  require_release(identical(file_md5(source), recorded$md5),
    "Frozen definitions do not reproduce the original checksum; restore the original RDS before release.")
  list(definitions = readRDS(source), source = source, recovered = recovered,
       evidence = evidence, md5 = recorded$md5)
}

resolve_validation_inputs <- function(vinputs, rows, root, proc, definitions) {
  original_root <- project_path(proc$output_dir, proc$project_dir)
  prefix <- paste0(sub("/+$", "", original_root), "/")
  resolved <- vinputs$path[rows]
  actions <- rep("Original source path", length(rows))
  for (i in seq_along(rows)) {
    original <- vinputs$path[rows[i]]
    if (startsWith(original, prefix)) {
      relative <- substring(original, nchar(prefix) + 1L)
      if (identical(relative, "logs/scoring_definitions.rds")) {
        resolved[i] <- definitions$source
        actions[i] <- if (definitions$recovered) "Byte-identical recovery from archived per-array scoring details" else "Frozen definitions in current processing folder"
      } else {
        resolved[i] <- owned_file(root, relative)
        actions[i] <- "Verified relative path in current processing folder"
      }
    }
  }
  require_release(all(file.exists(resolved)), "A recorded validation input is missing.")
  observed <- file_md5(resolved)
  require_release(identical(observed, vinputs$md5[rows]),
                  "A validation input changed after Step 03; validate again before release.")
  data.frame(recorded_path = vinputs$path[rows], resolved_path = resolved,
    recorded_md5 = vinputs$md5[rows], observed_md5 = observed,
    resolution = actions, passed = TRUE, stringsAsFactors = FALSE)
}

read_overlay_context <- function(context, config) {
  qa <- config$overlay_qa_dir
  if (is.null(qa)) {
    base <- file.path(config$project_dir, paste0("spatial_overlay_QA_v", context$version))
    candidates <- c(list.dirs(base, recursive = FALSE, full.names = TRUE),
                    file.path(context$root, "verification", "spatial_overlays"))
    candidates <- candidates[file.exists(file.path(candidates, "qa_run_summary.csv"))]
    valid <- lapply(candidates, function(p) {
      x <- tryCatch(fread(file.path(p, "qa_run_summary.csv")), error = function(e) NULL)
      if (is.null(x) || nrow(x) != 1L || !all(c("processing_run_id", "run_utc") %in% names(x)) ||
          !identical(x$processing_run_id, context$proc$run_id) ||
          !file.exists(file.path(p, "RUN_STATUS.txt")) ||
          !any(startsWith(readLines(file.path(p, "RUN_STATUS.txt")), "COMPLETE:"))) return(NULL)
      data.frame(path = p, run_utc = x$run_utc)
    })
    valid <- rbindlist(valid)
    require_release(nrow(valid) > 0L, "Run Step 04 first, or supply --overlay-dir=PATH to one completed QA run.")
    valid <- valid[order(run_utc, path, decreasing = TRUE)]
    qa <- valid$path[1]
  }
  qa <- normalizePath(project_path(qa, config$project_dir), winslash = "/", mustWork = TRUE)
  message("Selected Step 04 QA run: ", qa)
  summary <- fread(owned_file(qa, "qa_run_summary.csv"))
  inventory <- fread(owned_file(qa, "automatic_coordinate_checks_all_arrays.csv"))
  selection <- fread(owned_file(qa, "representative_selection.csv"))
  review <- fread(owned_file(qa, "overlay_review.csv"))
  inputs <- fread(owned_file(qa, "overlay_input_manifest.csv"))
  needed <- c("subarray_id", "count_id", "n_spots", "automatic_geometry_ok",
    "reflected_relative_to_grid", "source_spot_set_matches", "exported_coordinates_match_source",
    "max_source_coordinate_difference", "n_nonfinite_pixel_coordinates", "n_spots_outside_image", "affine_max_residual_px")
  require_release(all(needed %in% names(inventory)) && nrow(inventory) == nrow(context$good) &&
    !anyDuplicated(inventory$subarray_id) && setequal(inventory$subarray_id, context$good$subarray_id),
    "Overlay coordinate audit does not cover exactly the successful processing cohort.")
  require_release(nrow(summary) == 1L && identical(summary$processing_run_id, context$proc$run_id) &&
    canonical_release_version(summary$release_version) == context$version &&
    summary$n_arrays_checked == nrow(inventory) &&
    any(startsWith(readLines(owned_file(qa, "RUN_STATUS.txt")), "COMPLETE:")),
    "Overlay QA has a different run/version or is incomplete.")
  i <- match(inventory$subarray_id, context$good$subarray_id)
  require_release(all(inventory$n_spots == context$good$n_spots[i]) &&
    all(as.character(inventory$count_id) == as.character(context$good$count_id[i])),
    "Overlay array/specimen/spot counts disagree with processing.")
  require_release(all(c("subarray_id", "count_id", "figure_file") %in% names(selection)) &&
    nrow(selection) > 0L && !anyDuplicated(selection$subarray_id) &&
    all(selection$subarray_id %in% inventory$subarray_id) &&
    all(selection$figure_file == paste0("figures/", selection$subarray_id, "_overlay.png")),
    "Invalid representative selection or figure paths.")
  require_release(all(c("subarray_id", "image_alignment_review", "rotation_flip_review",
    "manual_annotation_review") %in% names(review)) &&
    !anyDuplicated(review$subarray_id) && setequal(review$subarray_id, selection$subarray_id),
    "Review rows must match the representative selection exactly.")
  review <- review[match(selection$subarray_id, review$subarray_id)]
  for (col in c("image_alignment_review", "rotation_flip_review", "manual_annotation_review")) {
    review[[col]] <- tolower(trimws(review[[col]]))
  }
  require_release(all(review$image_alignment_review %in% c("pending", "pass", "fail", "uncertain")) &&
    all(review$rotation_flip_review %in% c("pending", "pass", "fail", "uncertain")) &&
    all(review$manual_annotation_review %in% c("pending", "pass", "fail", "uncertain", "not_available")),
    "Review decisions must be pass/fail/uncertain/pending (or not_available for manual annotation).")
  require_release(all(c("role", "path", "md5") %in% names(inputs)), "Invalid overlay input manifest.")
  # Tie every checked metadata file to the exact current processing payload.
  spot_inputs <- inputs[role %in% c("released_spot_metadata", "processed_spot_metadata")]
  expected_names <- basename(context$good$spot_table_rds_file)
  require_release(nrow(spot_inputs) == nrow(context$good) && !anyDuplicated(basename(spot_inputs$path)) &&
    setequal(basename(spot_inputs$path), expected_names), "Overlay input manifest has a different spot-table inventory.")
  observed <- file_md5(vapply(context$good$spot_table_rds_file, function(p) owned_file(context$root, p), character(1)))
  require_release(identical(observed, spot_inputs$md5[match(expected_names, basename(spot_inputs$path))]),
                  "Spot metadata changed after overlay QA; rerun Step 04.")
  script <- inputs[role == "overlay_script"]
  require_release(nrow(script) == 1L && identical(file_md5(owned_file(qa, "overlay_script.R")), script$md5),
                  "The archived overlay script differs from the script used for QA.")
  image_inputs <- inputs[role == "source_images"]
  expected_images <- context$inputs[input %in% paste0("images_", unique(context$good$count_id))]
  require_release(nrow(image_inputs) == nrow(expected_images) &&
    !anyDuplicated(basename(image_inputs$path)) &&
    setequal(basename(image_inputs$path), basename(expected_images$path)) &&
    identical(image_inputs$md5[match(basename(expected_images$path), basename(image_inputs$path))], expected_images$md5) &&
    identical(file_md5(expected_images$path), expected_images$md5),
    "Overlay source images differ from the frozen processing inputs.")
  lookup <- inputs[role == "source_array_lookup"]
  ids <- context$inputs[input == "ids"]
  require_release(nrow(lookup) == 1L && identical(lookup$md5, ids$md5), "Overlay source ID lookup differs from processing.")
  if (file.exists(file.path(qa, "qa_output_manifest.csv"))) {
    outputs <- fread(file.path(qa, "qa_output_manifest.csv"))
    require_release(all(c("relative_path", "md5", "file_size_bytes") %in% names(outputs)) &&
      !anyDuplicated(outputs$relative_path), "Invalid QA output manifest.")
    protected <- c("RUN_STATUS.txt", "qa_run_summary.csv", "automatic_coordinate_checks_all_arrays.csv",
      "representative_selection.csv", "overlay_input_manifest.csv", "overlay_script.R", "sessionInfo.txt",
      "README_overlay_check.txt", "overlay_contact_sheet.png", selection$figure_file)
    require_release(all(protected %in% outputs$relative_path) &&
      !any(outputs$relative_path %in% c("overlay_review.csv", "qa_output_manifest.csv")),
      "The QA output manifest does not protect the required generated outputs.")
    paths <- vapply(outputs$relative_path, function(p) owned_file(qa, p), character(1))
    require_release(identical(file_md5(paths), outputs$md5) &&
      all(file.info(paths)$size == outputs$file_size_bytes), "A protected QA output changed; rerun Step 04.")
  }
  figures <- c(selection$figure_file, "overlay_contact_sheet.png")
  png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  for (figure in figures) require_release(
    identical(readBin(owned_file(qa, figure), "raw", n = 8L), png_signature), paste("Missing/invalid PNG:", figure))
  coordinate_ok <- all(inventory$automatic_geometry_ok %in% TRUE) &&
    all(inventory$source_spot_set_matches %in% TRUE) && all(inventory$exported_coordinates_match_source %in% TRUE) &&
    all(is.finite(inventory$max_source_coordinate_difference)) && all(inventory$max_source_coordinate_difference <= 1e-6) &&
    all(inventory$n_nonfinite_pixel_coordinates == 0L) && all(inventory$n_spots_outside_image == 0L) &&
    all(is.finite(inventory$affine_max_residual_px)) && all(inventory$affine_max_residual_px <= 1)
  coverage_ok <- nrow(selection) >= 8L && all(inventory$subarray_id[inventory$reflected_relative_to_grid %in% TRUE] %in% selection$subarray_id)
  visual_pass <- review$image_alignment_review %in% "pass" & review$rotation_flip_review %in% "pass"
  blockers <- c(if (!coordinate_ok) "Overlay coordinate checks have unresolved flags.",
    if (!coverage_ok) "Representative QA requires at least 8 arrays and all reflected source cases.",
    if (any(review$manual_annotation_review %in% "fail")) "Manual annotation review records an unresolved failure.",
    if (!all(visual_pass)) paste0(sum(!visual_pass), " representative arrays lack passing alignment/orientation decisions."))
  names <- c("RUN_STATUS.txt", "qa_run_summary.csv", "automatic_coordinate_checks_all_arrays.csv",
    "representative_selection.csv", "overlay_review.csv", "overlay_input_manifest.csv", "overlay_script.R",
    "sessionInfo.txt", "README_overlay_check.txt", "qa_output_manifest.csv", figures)
  files <- names[file.exists(file.path(qa, names))]
  list(directory = qa, summary = summary, inventory = inventory, selection = selection, review = review,
    ready = !length(blockers), blockers = blockers, coordinate_ok = coordinate_ok,
    coverage_ok = coverage_ok, n_visual_pass = sum(visual_pass), visual_pass = visual_pass,
    visual_status = if (all(visual_pass)) "passed" else "incomplete_or_flagged",
    manual_status = if (any(review$manual_annotation_review %in% "fail")) "flagged" else
      if (all(review$manual_annotation_review %in% c("pass", "not_available"))) "passed_where_available" else "not_completed_or_uncertain",
    files = unique(files), protected_output_manifest = file.exists(file.path(qa, "qa_output_manifest.csv")))
}

numeric_equal <- function(x, y, tolerance = 1e-10) {
  if (!is.numeric(x) || length(x) != length(y)) return(FALSE)
  x <- as.numeric(x); y <- as.numeric(y)
  identical(is.na(x), is.na(y)) &&
    all(is.finite(x[!is.na(x)])) &&
    all(abs(x[!is.na(y)] - y[!is.na(y)]) <= tolerance)
}

reference_z <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  values <- x[ok]
  center <- median(values)
  spread <- stats::mad(values, center = center, constant = 1.4826)
  if (is.finite(spread) && spread > 0) {
    out[ok] <- (values - center) / spread
  } else {
    spread <- stats::sd(values)
    out[ok] <- if (is.finite(spread) && spread > 0) {
      (values - mean(values)) / spread
    } else rep(0, length(values))
  }
  out
}

stamp_release <- function(x, context) {
  x <- as.data.frame(x)
  x$release_version <- context$version
  x$release_date <- context$date
  x$processing_run_id <- context$proc$run_id
  x$validation_run_id <- context$val$validation_run_id
  x
}

read_release_context <- function(config, version, script_file) {
  release_script_md5 <- file_md5(script_file)
  root <- resolve_processing_root(config, version)
  val_path <- project_path(config$validation_results_dir, config$project_dir)
  if (!file.exists(file.path(val_path, "validation_config.rds"))) val_path <- file.path(root, "abundance_validation")
  val_dir <- normalizePath(val_path, winslash = "/", mustWork = TRUE)
  require_release(startsWith(val_dir, paste0(root, "/")),
                  "Validation must belong to this processing release.")
  proc <- readRDS(file.path(root, "logs", "processing_config.rds"))
  val <- readRDS(file.path(val_dir, "validation_config.rds"))
  pr <- fread(file.path(root, "logs", "processing_run_summary.csv"))
  vr <- fread(file.path(val_dir, "validation_run_summary.csv"))
  require_release(
    nrow(pr) == 1L && all(c("run_status", "run_id", "full_cohort") %in% names(pr)) &&
      all(pr$run_status %in% "complete") && all(pr$full_cohort %in% TRUE) &&
      all(pr$run_id %in% proc$run_id),
    "A completed full-cohort Step 02 run is required."
  )
  require_release(
    nrow(vr) == 1L &&
      all(c("validation_status", "run_id", "validation_run_id") %in% names(vr)) &&
      all(vr$validation_status %in% "complete") &&
      all(vr$run_id %in% proc$run_id) &&
      all(vr$validation_run_id %in% val$validation_run_id) &&
      identical(val$run_id, proc$run_id),
    "A completed Step 03 run with matching identifiers is required."
  )
  status <- fread(file.path(root, "status", "subarray_processing_status_latest.csv"))
  qc <- fread(file.path(root, "cohort_qc_summary.csv"))
  coverage <- fread(file.path(root, "cohort_module_coverage.csv"))
  eligibility <- fread(file.path(val_dir, "subarray_eligibility_for_review.csv"))
  for (x in list(status, qc, coverage, eligibility)) {
    require_release("run_id" %in% names(x) && all(x$run_id %in% proc$run_id),
                    "Processing/validation tables mix run identifiers.")
  }
  require_release(all(eligibility$validation_run_id %in% val$validation_run_id),
                  "Eligibility records mix validation identifiers.")
  inputs <- fread(file.path(root, "logs", "processing_input_manifest.csv"))
  ids_entry <- inputs[input == "ids"]
  require_release(nrow(ids_entry) == 1L &&
                    identical(file_md5(ids_entry$path), ids_entry$md5),
                  "The source ID object differs from the processing input.")
  ids <- readRDS(ids_entry$path)
  require_release(
    nrow(ids) == 281L && !anyDuplicated(rownames(ids)) &&
      nrow(status) == 281L && !anyDuplicated(status$subarray_id) &&
      setequal(rownames(ids), status$subarray_id),
    "The source/status inventory must contain exactly 281 unique arrays."
  )
  good <- status[success %in% TRUE]
  bad <- status[!success %in% TRUE]
  require_release(
    nrow(good) == 280L && nrow(bad) == 1L &&
      identical(bad$subarray_id, "CN17_E2") &&
      !is.na(bad$error) &&
      grepl("No spots found for subarray: CN17_E2", bad$error, fixed = TRUE),
    "Unexpected processing failures or missing arrays remain."
  )
  require_release(
    nrow(qc) == 280L && !anyDuplicated(qc$subarray_id) &&
      setequal(qc$subarray_id, good$subarray_id) &&
      sum(qc$n_spots) == 270310L && uniqueN(qc$count_id) == 94L &&
      uniqueN(qc$subject_id) == 92L &&
      nrow(eligibility) == 281L && !anyDuplicated(eligibility$subarray_id) &&
      setequal(eligibility$subarray_id, status$subarray_id),
    "Cohort counts or eligibility records are inconsistent."
  )
  expected_subject <- as.character(qc$count_id)
  expected_subject[qc$count_id == 30] <- "53"
  expected_subject[qc$count_id == 66] <- "58"
  require_release(
    identical(as.character(qc$subject_id), expected_subject) &&
      all(qc$specimen_type == ifelse(qc$count_id %in% c(30, 66), "relapse", "primary")),
    "Patient/specimen linkage is inconsistent."
  )
  vinputs <- fread(file.path(val_dir, "validation_input_manifest.csv"))
  definitions <- read_frozen_definitions(root, good, vinputs)
  defs <- definitions$definitions
  modules <- names(defs$gene_sets)
  require_release(
    length(modules) == 19L && nrow(coverage) == 280L * length(modules) &&
      !anyDuplicated(coverage[, .(subarray_id, module)]) &&
      all(coverage$subarray_id %in% good$subarray_id) &&
      all(coverage$module %in% modules),
    "Module coverage is incomplete or duplicated."
  )

  # Verify archived code, allowing subsequent edits to working scripts.
  for (role in c("config", "runner", "helpers", "gene_sets")) {
    recorded <- inputs[input == role]
    archive <- file.path(root, "logs", "source_code", paste0(role, ".R"))
    require_release(nrow(recorded) == 1L &&
                      identical(file_md5(archive), recorded$md5),
                    paste("Processing source snapshot mismatch:", role))
  }
  require_release(all(vinputs$run_id %in% proc$run_id) &&
                    all(vinputs$validation_run_id %in% val$validation_run_id),
                  "Validation input manifest has inconsistent identifiers.")
  script_row <- which(basename(vinputs$path) == "03_validate_wang_marker_scores.R")
  config_row <- which(basename(vinputs$path) == "00_wang_config.R")
  require_release(length(script_row) == 1L,
                  "The validation input manifest must identify exactly one validation script.")
  validation_script_md5 <- file_md5(file.path(val_dir, "validation_script.R"))
  validation_script_verified <- identical(validation_script_md5, vinputs$md5[script_row])
  provenance_checks <- data.frame(
    check = "archived_validation_script_matches_recorded_input",
    passed = validation_script_verified,
    recorded_md5 = vinputs$md5[script_row],
    observed_md5 = validation_script_md5,
    stringsAsFactors = FALSE
  )
  data_rows <- setdiff(seq_len(nrow(vinputs)), c(script_row, config_row))
  input_resolution <- resolve_validation_inputs(vinputs, data_rows, root, proc, definitions)
  provenance_checks <- rbind(provenance_checks, data.frame(
    check = "frozen_scoring_definitions_match_original_checksum", passed = TRUE,
    recorded_md5 = definitions$md5, observed_md5 = file_md5(definitions$source)))

  summary <- fread(file.path(val_dir, "predictor_performance_summary.csv"))
  require_release(
    all(summary$run_id %in% proc$run_id) &&
      all(summary$validation_run_id %in% val$validation_run_id) &&
      setequal(summary$comparison_set, c("common_predictors", "available_predictor")) &&
      setequal(summary$prediction_type, c(
        "source_specimen_LOO", "source_full_fit", "fixed_module_annotation_agreement"
      )),
    "Validation summaries do not contain the required comparison/provenance groups."
  )
  list(
    root = root, val_dir = val_dir, proc = proc, val = val,
    status = status, good = good, qc = qc, coverage = coverage,
    eligibility = eligibility, ids = ids, defs = defs, inputs = inputs,
    definitions = definitions, input_resolution = input_resolution,
    version = version, date = format(Sys.time(), "%Y-%m-%d", tz = "America/New_York"),
    script_file = script_file, release_script_md5 = release_script_md5,
    validation_script_verified = validation_script_verified,
    provenance_checks = provenance_checks
  )
}

make_file_plan <- function(context) {
  rows <- list()
  add <- function(relative, role, source = NULL) {
    relative <- relative_release_path(relative)
    if (is.null(source)) source <- owned_file(context$root, relative)
    require_release(file.exists(source), paste("Missing release input:", source))
    rows[[length(rows) + 1L]] <<- data.frame(
      relative_path = relative, source = normalizePath(source, winslash = "/"),
      role = role, stringsAsFactors = FALSE
    )
  }
  for (i in seq_len(nrow(context$good))) {
    row <- context$good[i]
    sid <- row$subarray_id
    add(row$normalized_expression_file, "normalized_expression")
    add(row$spot_table_rds_file, "spot_rds")
    if (!is.na(row$spot_table_csv_file) && nzchar(row$spot_table_csv_file)) {
      add(row$spot_table_csv_file, "spot_csv")
    }
    if (!is.na(row$object_file) && nzchar(row$object_file)) {
      add(row$object_file, "seurat")
    }
    add(file.path("qc_summaries", paste0(sid, "_qc.csv")), "qc")
    add(file.path("module_coverage", paste0(sid, "_module_coverage.csv")), "coverage")
    add(file.path("module_coverage", paste0(sid, "_scoring_details.rds")), "scoring_details")
  }
  for (name in c(
    "cohort_qc_summary.csv", "cohort_module_coverage.csv",
    "status/subarray_processing_status.csv",
    "status/subarray_processing_status_latest.csv",
    "logs/processing_config.rds", "logs/processing_input_manifest.csv",
    "logs/processing_run_summary.csv", "logs/processing_session_info.txt",
    paste0("logs/source_code/", c("config", "runner", "helpers", "gene_sets"), ".R")
  )) add(name, "processing_provenance")
  add("logs/scoring_definitions.rds", "frozen_scoring_definitions", context$definitions$source)

  validation_files <- c(
    "annotated_spots.rds", "annotated_validation_spots.rds",
    "candidate_predictor_median_auc.pdf", "candidate_predictor_median_rho.pdf",
    "cohort_validation_summary.csv", "predictor_performance_by_subarray.csv",
    "predictor_performance_summary.csv", "subarray_eligibility_for_review.csv",
    "validation_config.rds", "validation_input_manifest.csv",
    "validation_run_summary.csv", "validation_script.R", "validation_session_info.txt"
  )
  for (name in validation_files) {
    add(file.path("abundance_validation", name), "validation",
        owned_file(context$val_dir, name))
  }
  for (name in context$overlay$files) {
    add(file.path("verification", "spatial_overlays", name), "overlay_qa",
        owned_file(context$overlay$directory, name))
  }
  plan <- rbindlist(rows)
  require_release(!anyDuplicated(plan$relative_path), "Duplicate release file destinations.")
  plan$source_md5 <- file_md5(plan$source)
  plan
}

read_array_metadata <- function(context, row) {
  m <- as.data.frame(readRDS(owned_file(context$root, row$spot_table_rds_file)),
                     check.names = FALSE)
  sid <- row$subarray_id
  require_release("spot_id" %in% names(m) &&
                    !anyNA(m$spot_id) && !anyDuplicated(m$spot_id) &&
                    all(startsWith(m$spot_id, paste0(sid, "."))),
                  paste("Invalid spot identifiers:", sid))
  repaired <- !"subarray_id" %in% names(m)
  if (repaired) m$subarray_id <- sid
  required <- c(
    "run_id", "count_id", "subject_id", "specimen_type", "baseline_primary_specimen",
    "x_um", "y_um", "nCount_RNA", "nFeature_RNA", "qc_pass",
    "has_manual_annotation", "manual_annotated_pixels",
    "manual_annotation_supported", "manual_annotation_min_pixels",
    "counts_branch", "module_score_version"
  )
  require_release(all(required %in% names(m)) && nrow(m) == row$n_spots &&
                    all(m$subarray_id %in% sid) &&
                    all(m$run_id %in% context$proc$run_id) &&
                    all(m$count_id %in% row$count_id),
                  paste("Wrong metadata schema/array/run:", sid))
  require_release(all(vapply(m, function(x) {
    is.atomic(x) && is.null(dim(x)) && length(x) == nrow(m)
  }, logical(1))) && all(is.finite(m$x_um)) && all(is.finite(m$y_um)),
  paste("Non-scalar metadata or invalid coordinates:", sid))
  q <- context$qc[match(sid, context$qc$subarray_id)]
  require_release(all(m$subject_id %in% q$subject_id) &&
                    all(m$specimen_type %in% q$specimen_type) &&
                    all(m$baseline_primary_specimen %in% q$baseline_primary_specimen) &&
                    all(m$counts_branch %in% "countsNonCorrected") &&
                    all(is.finite(m$nCount_RNA)) && all(m$nCount_RNA > 0),
                  paste("Incorrect patient linkage/count provenance:", sid))
  threshold <- context$proc$min_annotation_biological_pixels
  supported <- is.finite(m$manual_annotated_pixels) & m$manual_annotated_pixels > threshold
  require_release(all(m$manual_annotation_min_pixels %in% threshold) &&
                    identical(unname(as.logical(m$manual_annotation_supported)), unname(supported)),
                  paste("Incorrect annotation support:", sid))
  fractions <- grep("^manual_fraction_", names(m), value = TRUE)
  zero <- is.finite(m$manual_annotated_pixels) & m$manual_annotated_pixels == 0
  if (length(fractions) && any(zero)) {
    require_release(all(is.na(as.matrix(m[zero, fractions, drop = FALSE]))),
                    paste("Unannotated spots have nonmissing fractions:", sid))
  }
  list(table = m, repaired = repaired)
}

schema_rows <- function(x, table) {
  require_release(all(vapply(x, function(z) {
    is.atomic(z) && is.null(dim(z)) && length(z) == nrow(x)
  }, logical(1))), paste("Invalid metadata columns:", table))
  data.frame(
    variable = names(x),
    data_type = vapply(x, function(z) paste(class(z), collapse = "/"), character(1)),
    table = table, stringsAsFactors = FALSE
  )
}

inspect_release_arrays <- function(context) {
  checks <- schemas <- mappings <- vector("list", nrow(context$good))
  fixed_sets <- lapply(context$defs$gene_sets, function(g) unique(as.character(g)))
  for (i in seq_len(nrow(context$good))) {
    row <- context$good[i]
    sid <- row$subarray_id
    metadata <- read_array_metadata(context, row)
    m <- metadata$table
    expr <- readRDS(owned_file(context$root, row$normalized_expression_file))
    details <- readRDS(owned_file(context$root,
      file.path("module_coverage", paste0(sid, "_scoring_details.rds"))))
    info <- details$module_scoring
    require_release(
      inherits(expr, "sparseMatrix") &&
        identical(colnames(expr), m$spot_id) && ncol(expr) == nrow(m) &&
        !anyDuplicated(rownames(expr)) && all(is.finite(expr@x)) && all(expr@x >= 0) &&
        identical(details$run_id, context$proc$run_id) &&
        identical(info$gene_sets, fixed_sets) &&
        identical(info$alias_groups, context$defs$alias_groups) &&
        identical(info$missing_policy, "NA_if_any_requested_gene_unavailable") &&
        identical(details$feature_map$assay_gene, rownames(expr)) &&
        all(m$module_score_version %in% info$version),
      paste("Expression/scoring definition mismatch:", sid)
    )
    require_release(all(m$nFeature_RNA == Matrix::colSums(expr > 0)),
                    paste("Detected-gene counts differ from expression:", sid))
    # Check raw-count log normalization and recovered UMI arithmetic.
    reconstructed <- expr
    reconstructed@x <- expm1(reconstructed@x)
    reconstructed <- reconstructed %*% Matrix::Diagonal(x = m$nCount_RNA / 10000)
    require_release(
      all(abs(Matrix::colSums(reconstructed) - m$nCount_RNA) < 1e-6) &&
        all(abs(reconstructed@x - round(reconstructed@x)) < 1e-6),
      paste("Raw-UMI normalization check failed:", sid)
    )
    rm(reconstructed)
    require_release(
      nrow(info$gene_mapping) == sum(lengths(fixed_sets)) &&
        !anyDuplicated(info$gene_mapping[, c("module", "requested_gene")]),
      paste("Incomplete gene mapping:", sid)
    )
    for (module in names(fixed_sets)) {
      map <- info$gene_mapping[info$gene_mapping$module == module, , drop = FALSE]
      require_release(identical(map$requested_gene, fixed_sets[[module]]),
                      paste("Wrong marker definition:", sid, module))
      complete <- !anyNA(map$assay_gene)
      require_release(!anyDuplicated(map$assay_gene[!is.na(map$assay_gene)]) &&
                        all(map$assay_gene[!is.na(map$assay_gene)] %in% rownames(expr)),
                      paste("Invalid/duplicated gene mapping:", sid, module))
      expected <- if (complete) {
        as.numeric(Matrix::colMeans(expr[map$assay_gene, , drop = FALSE]))
      } else rep(NA_real_, nrow(m))
      col <- paste0("score_", module)
      require_release(all(c(col, paste0(col, "_z"), paste0(col, "_pos")) %in% names(m)) &&
                        numeric_equal(m[[col]], expected) &&
                        numeric_equal(m[[paste0(col, "_z")]], reference_z(expected)) &&
                        numeric_equal(m[[paste0(col, "_pos")]], pmax(reference_z(expected), 0)),
                      paste("Module mean/scaling/missingness mismatch:", sid, module))
      selected <- which(context$coverage$subarray_id == sid &
                          context$coverage$module == module)
      cv <- context$coverage[selected]
      present <- !is.na(map$assay_gene)
      n_detected <- if (any(present)) {
        sum(Matrix::rowSums(expr[map$assay_gene[present], , drop = FALSE] > 0) > 0)
      } else 0L
      require_release(nrow(cv) == 1L && cv$n_requested == length(fixed_sets[[module]]) &&
                        cv$score_denominator == length(fixed_sets[[module]]) &&
                        cv$n_available == sum(!is.na(map$assay_gene)) &&
                        cv$n_detected == n_detected &&
                        abs(cv$fraction_available - mean(present)) < 1e-10 &&
                        abs(cv$fraction_detected - n_detected / nrow(map)) < 1e-10 &&
                        cv$score_eligible == complete,
                      paste("Module coverage mismatch:", sid, module))
    }
    mappings[[i]] <- data.frame(
      subarray_id = sid, count_id = as.character(row$count_id),
      info$gene_mapping, stringsAsFactors = FALSE
    )
    schemas[[i]] <- schema_rows(m, "spot_tables")
    checks[[i]] <- data.frame(
      subarray_id = sid, n_spots = nrow(m),
      expression_metadata_aligned = TRUE, normalization_verified = TRUE,
      module_scores_verified = TRUE, annotation_support_verified = TRUE,
      subarray_id_restored = metadata$repaired, stringsAsFactors = FALSE
    )
    rm(expr, details, m)
    if (i %% 40L == 0L) message("Verified ", i, "/", nrow(context$good), " arrays")
  }
  for (name in c("annotated_spots.rds", "annotated_validation_spots.rds")) {
    x <- as.data.frame(readRDS(file.path(context$val_dir, name)))
    require_release(all(x$run_id %in% context$proc$run_id) &&
                      all(x$validation_run_id %in% context$val$validation_run_id),
                    paste("Validation RDS has inconsistent identifiers:", name))
    schemas[[length(schemas) + 1L]] <- schema_rows(x, name)
  }
  list(checks = rbindlist(checks), schema = rbindlist(schemas),
       mapping = rbindlist(mappings))
}

build_cohort_manifest <- function(context) {
  ids <- context$ids
  x <- data.frame(
    subarray_id = rownames(ids), count_id = as.character(ids$id),
    parent_slide = as.character(ids$slide), subslide = as.character(ids$subslide),
    sample_label = as.character(ids$ID),
    source_has_manual_annotation = as.logical(ids$hasAnnot),
    stringsAsFactors = FALSE
  )
  s <- context$status[match(x$subarray_id, context$status$subarray_id)]
  x$processing_success <- s$success
  x$processing_timestamp <- as.character(s$timestamp)
  x$processing_error <- s$error
  for (col in c("normalized_expression_file", "spot_table_rds_file",
                "spot_table_csv_file", "object_file")) x[[col]] <- s[[col]]
  names(x)[names(x) == "object_file"] <- "processed_object_file"
  q <- context$qc[match(x$subarray_id, context$qc$subarray_id)]
  for (col in setdiff(names(q), c("run_id", "count_id", "subarray_id"))) x[[col]] <- q[[col]]
  e <- context$eligibility[match(x$subarray_id, context$eligibility$subarray_id)]
  for (col in c("subject_id", "specimen_type", "baseline_primary_specimen",
                "eligible_spatial_main", "eligible_baseline_primary", "eligible_primary",
                "eligible_tumor", "eligible_lymphocyte", "eligible_tls", "review_status")) {
    require_release(col %in% names(e), paste("Missing eligibility field:", col))
    x[[col]] <- e[[col]]
  }
  x$inclusion_status <- ifelse(x$processing_success %in% TRUE,
                              "included", "excluded_no_retained_spots")
  x$scoring_details_file <- ifelse(x$processing_success,
    file.path("module_coverage", paste0(x$subarray_id, "_scoring_details.rds")), NA_character_)
  x <- x[order(as.numeric(x$count_id), x$subarray_id), ]
  stamp_release(x, context)
}

describe_variable <- function(variable, threshold, positive_threshold) {
  description <- switch(variable,
    spot_id = "Unique spot identifier; equals its expression-matrix column name.",
    count_id = "Source specimen identifier; 94 specimens correspond to 92 patients.",
    subarray_id = "Array identifier; any release-time restoration is recorded in metadata_repairs.csv.",
    subject_id = "Independent patient ID, with relapse 30 linked to 53 and relapse 66 linked to 58.",
    specimen_type = "Primary or relapse specimen; relapse sections are not independent baseline patients.",
    baseline_primary_specimen = "TRUE for a primary specimen, independent of spatial/follow-up eligibility.",
    run_id = "Identifier of the frozen Step 02 processing run.",
    validation_run_id = "Identifier of the completed Step 03 validation run.",
    subarray_id_reconstructed = "Historical flag: Step 03 reconstructed this identifier from its input filename.",
    counts_branch = "countsNonCorrected denotes source raw UMI input.",
    module_score_version = "Version of the frozen marker-mean calculation.",
    nCount_RNA = "Source raw UMI total for the spot.",
    nFeature_RNA = "Number of genes with positive expression in the spot.",
    x_um = "Centered within-array x coordinate; 150/sqrt(2) microns per original lattice unit.",
    y_um = "Centered within-array y coordinate; not a common registered frame across sections.",
    pixel_x = "Source image x coordinate (rightwards); see image_provenance.csv for technical overlay review.",
    pixel_y = "Source image y coordinate (downwards); technical image review does not approve pathology masks or boundaries.",
    manual_annotated_pixels = "Number of biologically annotated pixels; excludes Nothing, holes and artefacts.",
    manual_has_annotated_pixels = "Whether the spot has any biological annotation pixels.",
    manual_annotation_supported = paste0("Whether biological pixels are strictly greater than ", threshold, "."),
    manual_annotation_min_pixels = "Recorded annotation-support threshold; the comparison is strict greater-than.",
    eligible_primary = "Legacy compatibility name for the spatial analysis tier; not a primary-specimen flag.",
    eligible_spatial_main = "Array meets the main spatial/QC criteria, irrespective of primary/relapse status.",
    eligible_baseline_primary = "Array meets spatial criteria and is a primary specimen.",
    qc_pass = "Recorded UMI/artifact QC flag; source-retained spots are preserved in the release.",
    linked_specimen_in_training = "For source LOO predictions, a linked specimen from patient 53 or 58 may remain in training.",
    prediction_type = "Distinguishes source full fit, source specimen LOO, and fixed-module annotation agreement.",
    comparison_set = "Common predictors use matched supported spots; available-predictor results describe coverage.",
    n_spots = "Number of retained spots, or actual metric input pairs in a validation performance table.",
    n_target_eligible = "QC-passing supported spots with a usable manual target before predictor availability.",
    n_positive = "Positive manual labels among the actual evaluated prediction/target pairs.",
    n_negative = "Negative manual labels among the actual evaluated prediction/target pairs.",
    rho = "Within-array Spearman correlation; undefined for constant targets or scores.",
    auc = "Within-array AUROC with prespecified direction; a constant finite score is 0.5 when both classes suffice.",
    NULL
  )
  if (is.null(description)) {
    description <- if (grepl("^hist_loo_", variable)) {
      "Author-supplied expression-model prediction held out by source specimen; not jointly by linked patient."
    } else if (grepl("^hist_score_", variable)) {
      "Expression-derived full-data histological prediction; annotation consistency includes training data."
    } else if (grepl("^hist_prediction_", variable)) {
      "Provenance of the source expression-derived histological predictions."
    } else if (grepl("^hist_", variable)) {
      "Summary of full-data expression-derived histological scores."
    } else if (grepl("^manual_fraction_", variable)) {
      "Fraction among biological annotation pixels; NA when no biological pixels exist. Validation additionally requires the support threshold."
    } else if (grepl("^manual_pixels_", variable)) {
      "Original annotation-pixel count for the named category."
    } else if (grepl("^manual_.*_positive$", variable)) {
      paste0("Supported manual target greater than ", positive_threshold, "; unsupported labels remain NA.")
    } else if (variable %in% c("manual_tumor", "manual_lymphocyte", "manual_tls")) {
      "Source-defined broad tissue target; unsupported annotations are NA, not negative labels."
    } else if (grepl("^score_.*_z$", variable)) {
      "Within-array robust z-score with SD fallback; finite constants become zero and missing values remain NA."
    } else if (grepl("^score_.*_pos$", variable)) {
      "Within-array z-score truncated at zero; missing values remain NA. Not a cell proportion."
    } else if (grepl("^score_", variable)) {
      "Mean log-normalized expression over the fixed, alias-resolved marker panel. An unavailable required gene makes the module NA."
    } else if (grepl("^clinical_.*_time$", variable)) {
      "Original source survival time in years; not administratively capped by preprocessing."
    } else if (grepl("^clinical_.*_status$", variable)) {
      "Original source survival event indicator; repeated across a specimen's spots."
    } else if (grepl("^clinical_", variable)) {
      "Source clinical value preserved without outcome modeling; repeated spots are not independent clinical observations."
    } else if (grepl("^qc_", variable)) {
      "Recorded spot quality-control flag."
    } else if (grepl("^percent\\.", variable)) {
      "Percentage of raw counts assigned to the indicated feature group."
    } else {
      "See the frozen processing/validation code and table-specific context."
    }
  }
  category <- if (grepl("^clinical_", variable)) "clinical" else if (
    grepl("^manual_", variable)) "manual_annotation" else if (
    grepl("^hist_", variable)) "expression_histology" else if (
    grepl("^score_", variable)) "marker_score" else "metadata_or_validation"
  data.frame(variable = variable, category = category, description = description,
             stringsAsFactors = FALSE)
}

build_metadata_dictionary <- function(schema, context) {
  types <- schema[, .(
    data_type = paste(sort(unique(data_type)), collapse = ";"),
    tables = paste(sort(unique(table)), collapse = ";")
  ), by = variable]
  descriptions <- rbindlist(lapply(types$variable, describe_variable,
    threshold = context$proc$min_annotation_biological_pixels,
    positive_threshold = context$val$manual_positive_threshold))
  merge(types, descriptions, by = "variable", sort = TRUE)
}

verify_release_manifest <- function(root) {
  manifest <- fread(file.path(root, "manifests", "file_manifest.csv"))
  require_release(!anyDuplicated(manifest$relative_path), "Duplicate manifest paths.")
  paths <- vapply(manifest$relative_path, function(p) owned_file(root, p), character(1))
  require_release(identical(file_md5(paths), manifest$md5) &&
                    all(file.info(paths)$size == manifest$file_size_bytes),
                  "Final release checksums or file sizes do not match.")
  actual <- list.files(root, recursive = TRUE, full.names = FALSE, all.files = TRUE)
  actual <- setdiff(actual, "manifests/file_manifest.csv")
  require_release(setequal(actual, manifest$relative_path),
                  "Final release inventory differs from its file manifest.")
  invisible(TRUE)
}

publish_release_directory <- function(staging, target) {
  # Stage on the same filesystem. Preserve the old folder, then roll it back
  # if exposing the verified replacement fails. Nothing is deleted.
  backup <- NULL
  if (file.exists(target)) {
    require_release(dir.exists(target), "The release target exists as a file.")
    backup <- tempfile(paste0(".", basename(target), "_backup_",
      format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_"), tmpdir = dirname(target))
    require_release(file.rename(target, backup), "Could not preserve existing release as a backup.")
  }
  if (!file.rename(staging, target)) {
    if (!is.null(backup)) require_release(file.rename(backup, target),
      paste("Could not restore the original release. It remains at:", backup))
    stop("Could not expose staged release; original release was restored.", call. = FALSE)
  }
  if (!is.null(backup)) message("Previous release retained at: ", backup)
  invisible(backup)
}

assemble_release <- function(context, plan, inspection, target) {
  require_release(context$validation_script_verified,
    paste("Archived validation code differs from its recorded input hash.",
          "Rerun saved Step 03 code in a fresh validation directory before release."))
  require_release(
    sum(plan$role == "normalized_expression") == nrow(context$good) &&
      sum(plan$role == "spot_rds") == nrow(context$good) &&
      sum(plan$role == "scoring_details") == nrow(context$good) &&
      nrow(inspection$checks) == nrow(context$good),
    "Release plan or array verification is incomplete."
  )
  require_release(context$overlay$ready, paste(context$overlay$blockers, collapse = " "))
  parent <- dirname(target)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(paste0(".", basename(target), "_staging_"), tmpdir = parent)
  require_release(dir.create(staging), "Cannot create release staging directory.")
  staging <- normalizePath(staging, winslash = "/", mustWork = TRUE)
  completed <- FALSE
  on.exit(if (!completed) message("Incomplete staging directory retained: ", staging), add = TRUE)
  repairs <- list()
  for (i in seq_len(nrow(plan))) {
    item <- plan[i]
    destination <- file.path(staging, item$relative_path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (item$role == "spot_csv") next
    require_release(file.copy(item$source, destination, overwrite = FALSE),
                    paste("Could not copy:", item$relative_path))
    require_release(identical(file_md5(destination), item$source_md5),
                    paste("Source changed or copying failed:", item$relative_path))
    if (item$role == "spot_rds") {
      sid <- sub("_spots\\.rds$", "", basename(item$relative_path))
      row <- context$good[match(sid, context$good$subarray_id)]
      staged_context <- context
      staged_context$root <- staging
      metadata <- read_array_metadata(staged_context, row)
      if (metadata$repaired) {
        saveRDS(metadata$table, destination, compress = TRUE)
        repairs[[length(repairs) + 1L]] <- data.frame(
          relative_path = item$relative_path, subarray_id = sid,
          action = "Restore subarray_id after spot-prefix/QC checks",
          source_md5 = item$source_md5, release_md5 = file_md5(destination),
          stringsAsFactors = FALSE
        )
      }
      csv <- row$spot_table_csv_file
      if (!is.na(csv) && nzchar(csv)) {
        # CSV exports are regenerated from the same authoritative metadata.
        fwrite(metadata$table, file.path(staging, relative_release_path(csv)))
      }
    }
  }
  repair_table <- if (length(repairs)) rbindlist(repairs) else data.table(
    relative_path = character(), subarray_id = character(), action = character(),
    source_md5 = character(), release_md5 = character()
  )
  dir.create(file.path(staging, "manifests"), showWarnings = FALSE)
  dir.create(file.path(staging, "documentation"), showWarnings = FALSE)
  dir.create(file.path(staging, "verification"), showWarnings = FALSE)
  cohort <- build_cohort_manifest(context)
  fwrite(cohort, file.path(staging, "manifests", "cohort_manifest.csv"))
  patient <- as.data.table(cohort)[, .(
    subject_id = subject_id[1], specimen_type = specimen_type[1],
    n_listed_sections = .N, n_processed_sections = sum(processing_success),
    n_retained_spots = sum(n_spots, na.rm = TRUE)
  ), by = count_id]
  fwrite(stamp_release(patient, context),
         file.path(staging, "manifests", "specimen_patient_inventory.csv"))
  fwrite(repair_table, file.path(staging, "manifests", "metadata_repairs.csv"))
  lineage <- copy(plan)
  lineage$release_action <- ifelse(lineage$role == "spot_csv",
    "CSV re-exported from verified RDS metadata", "Copied unchanged")
  lineage$release_action[lineage$relative_path %in% repair_table$relative_path] <-
    "Metadata-only identifier restoration"
  for (i in which(lineage$role == "spot_csv")) {
    row <- context$good[match(lineage$relative_path[i], context$good$spot_table_csv_file)]
    source_row <- match(row$spot_table_rds_file, plan$relative_path)
    lineage$source[i] <- plan$source[source_row]
    lineage$source_md5[i] <- plan$source_md5[source_row]
  }
  fwrite(lineage, file.path(staging, "manifests", "source_file_lineage.csv"))
  fwrite(context$input_resolution,
         file.path(staging, "verification", "validation_input_resolution.csv"))
  if (context$definitions$recovered) {
    lineage[role == "frozen_scoring_definitions", release_action := paste0(
      "Byte-identical recovery from ", context$definitions$evidence, "; original validation MD5 verified")]
    fwrite(lineage, file.path(staging, "manifests", "source_file_lineage.csv"))
  }
  review_receipt <- copy(context$overlay$review)
  review_receipt$technical_visual_pass <- context$overlay$visual_pass
  review_receipt$record_type <- "Technical image/spot review; no pathology boundary approval"
  review_receipt$processing_run_id <- context$proc$run_id
  review_receipt$assembled_date <- context$date
  fwrite(review_receipt, file.path(staging, "verification", "spatial_overlays", "review_status_at_release.csv"))

  genes <- rbindlist(lapply(names(context$defs$gene_sets), function(module) {
    markers <- unique(context$defs$gene_sets[[module]])
    data.frame(module = module, requested_gene = markers,
      score_denominator = length(markers), weight = 1 / length(markers))
  }))
  fwrite(genes, file.path(staging, "documentation", "module_gene_sets.csv"))
  fwrite(inspection$mapping,
         file.path(staging, "documentation", "module_gene_mapping_by_subarray.csv"))
  fwrite(build_metadata_dictionary(inspection$schema, context),
         file.path(staging, "documentation", "data_dictionary.csv"))
  image_provenance <- as.data.table(cohort)[, .(
    subarray_id, count_id, parent_slide, subslide,
    coordinate_frame = "array_local_centered",
    image_source = context$inputs$path[match(paste0("images_", count_id), context$inputs$input)],
    image_md5_at_processing = context$inputs$md5[match(paste0("images_", count_id), context$inputs$input)],
    annotation_source = context$inputs$path[match(paste0("annotation_", count_id), context$inputs$input)],
    registration_review = "Representative image/spot review recorded separately; full-resolution mask/boundary approval not established"
  )]
  coordinate_index <- match(image_provenance$subarray_id, context$overlay$inventory$subarray_id)
  visual_index <- match(image_provenance$subarray_id, context$overlay$review$subarray_id)
  image_provenance$coordinate_check_status <- ifelse(is.na(coordinate_index), "no_retained_spots",
    ifelse(context$overlay$inventory$automatic_geometry_ok[coordinate_index], "passed", "flagged"))
  image_provenance$image_alignment_review <- ifelse(is.na(visual_index), "not_selected_for_visual_review",
    context$overlay$review$image_alignment_review[visual_index])
  image_provenance$rotation_flip_review <- ifelse(is.na(visual_index), "not_selected_for_visual_review",
    context$overlay$review$rotation_flip_review[visual_index])
  image_provenance$spatial_qa_report <- "verification/spatial_overlays/review_status_at_release.csv"
  fwrite(image_provenance, file.path(staging, "documentation", "image_provenance.csv"))
  fwrite(inspection$checks,
         file.path(staging, "verification", "release_verification_by_subarray.csv"))
  checks <- data.frame(
    check = c("processing_complete", "validation_complete", "cohort_hierarchy",
              "source_validation_inputs_unchanged", "expression_metadata_alignment",
              "raw_UMI_normalization", "fixed_module_scores_and_missingness",
              "annotation_support", "successful_records_only", "spatial_coordinate_checks",
              "representative_overlay_coverage", "recorded_image_alignment_and_orientation_review"),
    passed = TRUE
  )
  fwrite(checks, file.path(staging, "verification", "release_verification_summary.csv"))
  fwrite(context$provenance_checks,
         file.path(staging, "verification", "release_provenance_checks.csv"))
  require_release(file.copy(context$script_file,
    file.path(staging, "logs", "release_script.R"), overwrite = FALSE),
    "Could not archive release code.")
  require_release(identical(
    file_md5(file.path(staging, "logs", "release_script.R")), context$release_script_md5
  ), "Release code changed during assembly; save the script and run again.")
  saveRDS(list(
    release_version = context$version, release_date = context$date,
    processing_run_id = context$proc$run_id,
    validation_run_id = context$val$validation_run_id,
    release_script_md5 = context$release_script_md5,
    source_processing_directory = context$root, target_directory = target,
    overlay_qa_source_directory = context$overlay$directory,
    overlay_review_md5 = file_md5(file.path(context$overlay$directory, "overlay_review.csv")),
    n_visual_reviewed = context$overlay$n_visual_pass,
    frozen_definitions_byte_identically_recovered = context$definitions$recovered
  ), file.path(staging, "logs", "release_config.rds"))
  capture.output(sessionInfo(), file = file.path(staging, "logs", "release_session_info.txt"))

  readme <- c(
    "# Wang TNBC ST processed dataset",
    "", paste0("Release version: ", context$version),
    "Source study: Wang et al., Nature Communications (2024). https://doi.org/10.1038/s41467-024-54145-w",
    paste0("Release date: ", context$date),
    paste0("Processing version: ", context$proc$release_version),
    paste0("Processing run: ", context$proc$run_id),
    paste0("Validation run: ", context$val$validation_run_id),
    "", "## Contents",
    paste0(nrow(cohort), " listed arrays; ", nrow(context$good), " processed arrays; ",
           format(sum(context$qc$n_spots), big.mark = ",", scientific = FALSE, trim = TRUE),
           " retained spots; ", uniqueN(cohort$count_id), " specimens from ",
           uniqueN(cohort$subject_id), " independent patients."),
    "Relapse specimen 30 belongs to patient 53; relapse 66 belongs to patient 58. CN17_E2 has no retained spots.",
    "The cohort manifest uses paths relative to this release. The specimen/patient inventory preserves all specimens.",
    "", "## Expression and scores",
    "Expression matrices are sparse genes by spots; columns exactly match spot_id in the corresponding metadata.",
    "Normalization is log(1 + 10000 * raw_count / total_raw_counts_for_the_spot), using countsNonCorrected.",
    "Measured zero-expression features are retained. Source-to-assay feature names and module mappings are saved per array in module_coverage/*_scoring_details.rds.",
    "Seurat-compatible feature names replace underscores and pipes with hyphens; the original names remain in the feature maps.",
    "Each marker score uses the same requested alias-resolved panel and fixed denominator. Any genuinely unavailable required feature makes the module NA.",
    "Within-array _z and _pos scores retain missingness. They describe local relative enrichment and are not comparable patient-level abundance or cell proportions.",
    "TLS_WANG_30 was developed in the source cohort; these mean-expression modules do not reproduce every published GSVA or weighted signature.",
    "", "## Annotation and validation",
    paste0("Raw manual fractions are defined for positive biological-pixel counts. Validation requires strictly more than ",
           context$proc$min_annotation_biological_pixels, " biological pixels; unsupported targets are NA."),
    "Main validation compares primary specimens using matched supported spots and matched evaluable sections. Each patient has one annotated primary section in this source.",
    "Available-predictor summaries separately report coverage. Consult evaluable counts, especially for rare TLS labels.",
    "hist_score_* predictions are expression-derived full-data fits. Their annotation agreement includes training data.",
    "hist_loo_* predictions are supplied by the authors and held out by source specimen, not jointly by linked patient.",
    "Patients 53 and 58 may have a linked specimen in source-LOO training. This release does not establish H&E-only performance or external clinical prediction.",
    "Clinical endpoints retain original times/statuses. TNBC3 has zero follow-up; prognosis requires a stated handling rule and patient-grouped analysis.",
    "Ten-year administrative censoring is an analysis decision, not a modification of the source exports.",
    "", "## Spatial interpretation",
    "x_um/y_um use the verified 150-micrometer lattice pitch and are independently centered within each array.",
    "They are not a shared registered coordinate frame for combining sections.",
    "Images are referenced in documentation/image_provenance.csv and are not embedded in Seurat objects or this release.",
    paste0("Step 04 coordinate consistency passed for ", nrow(context$overlay$inventory),
      " successful arrays. Recorded technical image/spot alignment and orientation review passed for ",
      context$overlay$n_visual_pass, " representative arrays, including all reflected source mappings."),
    "Consult verification/spatial_overlays/review_status_at_release.csv for the recorded decisions and any comments.",
    "Only selected arrays received representative visual review; coordinate checks cover all successful arrays.",
    "The QA run summary/figure captions preserve the initial generated pending status; review_status_at_release.csv records the actual final decisions.",
    "Per-spot manual tumor fractions are diagnostic overlays. Full-resolution annotation-mask registration and pathology boundary approval remain separate work.",
    "Partial coverage, gaps and tissue discontinuities require support-aware region selection in subsequent spatial analyses.",
    "Raw-versus-author-corrected feature sensitivity and repeat-section agreement of the final spatial contrast remain separate analyses.",
    "", "## Release integrity and provenance",
    "Only successful processing records and the completed validation run were selected; unrelated source-directory files were excluded.",
    "Processing and validation code/configuration snapshots remain as originally recorded. Provenance input paths identify the original local sources.",
    "Release CSVs were re-exported from verified RDS metadata. Any restored subarray_id fields are listed in manifests/metadata_repairs.csv; other values and expression matrices were preserved.",
    "subarray_id_reconstructed in validation tables describes historical recovery during Step 03, rather than the current release-table schema.",
    "The file manifest hashes every payload file except the file manifest itself. The final directory is exposed only after inventory and checksum verification.",
    "Frozen scoring definitions are included in logs/scoring_definitions.rds. Any byte-identical recovery and relocated validation input paths are recorded in verification/validation_input_resolution.csv and manifests/source_file_lineage.csv.",
    "Step 04 generates overlay QA before Step 05 assembles this folder. One chosen QA run is packaged under verification/spatial_overlays.",
    "Use verification/release_summary.csv and verification/release_verification_summary.csv to check completion."
  )
  writeLines(readme, file.path(staging, "documentation", "README_processed_data.md"))
  writeLines(c(
    paste0("# Changes in ", context$version), "",
    "Retained the verified raw-count normalization and source spots.",
    "Added fixed-denominator, alias-resolved module scores, availability/detection coverage, and preservation of unavailable scores as NA.",
    "Applied annotation-support checks and separated full-fit, source specimen-LOO, and fixed-module annotation agreement.",
    "Distinguished patients, primary/relapse specimens, and spatial versus baseline eligibility.",
    "Packaged frozen provenance and exact file selection with independent numerical and checksum verification.",
    "Added Step 04 representative H&E/spot QA, recorded technical visual decisions, input hashes and figures before final Step 05 assembly.",
    "Retained the Wang_TNBC_v1.1.0 directory structure and included frozen scoring definitions for reproducible reassembly."
  ), file.path(staging, "documentation", "CHANGELOG.md"))
  fwrite(stamp_release(data.frame(
    release_ready = TRUE, n_listed_arrays = nrow(cohort),
    n_processed_arrays = nrow(context$good), n_spots = sum(context$qc$n_spots),
    n_specimens = uniqueN(cohort$count_id), n_patients = uniqueN(cohort$subject_id),
    n_metadata_repairs = nrow(repair_table),
    scope = "Corrected preprocessing, annotation validation and representative technical image/spot overlay QA",
    n_spatial_arrays_checked = nrow(context$overlay$inventory),
    n_representative_arrays_reviewed = context$overlay$n_visual_pass,
    spatial_coordinate_checks_passed = context$overlay$coordinate_ok,
    spatial_visual_review_status = context$overlay$visual_status,
    manual_annotation_visual_review_status = context$overlay$manual_status,
    pathology_boundary_approved = FALSE,
    cross_array_registration_validated = FALSE,
    spatial_qa_run_utc = context$overlay$summary$run_utc
  ), context), file.path(staging, "verification", "release_summary.csv"))

  paths <- list.files(staging, full.names = TRUE, recursive = TRUE, all.files = TRUE)
  relative <- substring(paths, nchar(staging) + 2L)
  manifest <- data.frame(
    relative_path = relative, file_name = basename(paths),
    file_size_bytes = file.info(paths)$size, md5 = file_md5(paths),
    release_version = context$version, release_date = context$date,
    processing_run_id = context$proc$run_id,
    validation_run_id = context$val$validation_run_id, stringsAsFactors = FALSE
  )
  fwrite(manifest, file.path(staging, "manifests", "file_manifest.csv"))
  verify_release_manifest(staging)
  require_release(identical(file_md5(plan$source), plan$source_md5),
                  "An input changed during assembly; incomplete staging is retained.")
  backup <- publish_release_directory(staging, target)
  completed <- TRUE
  message("Verified release ready: ", target)
  invisible(list(release_dir = target, backup_dir = backup, summary = cohort, file_manifest = manifest))
}

run_prepare_release <- function(config, dry_run = FALSE, release_version = NULL) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(config$project_dir)
  version <- canonical_release_version(
    if (is.null(release_version)) config$release_version else release_version
  )
  target <- if (!is.null(config$release_output_dir)) config$release_output_dir else
    file.path(project_path(config$working_output_dir, config$project_dir), paste0("Wang_TNBC_v", version))
  if (!grepl("^/", target)) target <- file.path(getwd(), target)
  if (dir.exists(target)) target <- normalizePath(target, winslash = "/", mustWork = TRUE)
  script_file <- normalizePath("05_prepare_wang_release.R", mustWork = TRUE)
  context <- read_release_context(config, version, script_file)
  require_release(target == context$root || !startsWith(context$root, paste0(sub("/+$", "", target), "/")),
    "The target cannot contain the processing input directory.")
  context$overlay <- read_overlay_context(context, config)
  if (!dry_run) {
    require_release(context$validation_script_verified,
      paste("Archived Step 03 code differs from its recorded input hash.",
            "Rerun saved Step 03 code in a fresh validation directory before release."))
  }
  message("Checking completed runs and ", nrow(context$good), " arrays for release ", version)
  plan <- make_file_plan(context)
  inspection <- inspect_release_arrays(context)
  ready <- context$validation_script_verified && context$overlay$ready
  if (dry_run) {
    message("Numerical/data checks passed; release_ready: ", ready)
    if (!context$validation_script_verified) {
      message("BLOCKER: archived Step 03 code differs from its recorded input hash.",
              " Rerun saved Step 03 code in a fresh validation directory before release.")
    }
    for (blocker in context$overlay$blockers) message("BLOCKER: ", blocker)
    if (context$definitions$recovered) message("Frozen scoring definitions recovered byte-identically; original checksum verified.")
    message("No release files written. Target: ", target)
    return(invisible(list(
      release_dir = target, release_version = version,
      release_ready = ready,
      n_arrays = nrow(context$good), n_spots = sum(context$qc$n_spots),
      metadata_repairs_needed = sum(inspection$checks$subarray_id_restored),
      checks = inspection$checks, plan = plan,
      provenance_checks = context$provenance_checks,
      overlay_directory = context$overlay$directory,
      n_visual_pass = context$overlay$n_visual_pass,
      blockers = context$overlay$blockers
    )))
  }
  require_release(ready, paste(context$overlay$blockers, collapse = " "))
  assemble_release(context, plan, inspection, target)
}

if (sys.nframe() == 0L) {
  raw_args <- commandArgs(TRUE)
  if (any(raw_args %in% c("--help", "-h"))) {
    cat("Usage: Rscript 05_prepare_wang_release.R [--dry-run] [--overlay-dir=QA_RUN]\n",
        "       [--processing-dir=PATH] [--validation-dir=PATH] [--release-dir=PATH]\n",
        "Keeps outputs_wang_v1.1/Wang_TNBC_v1.1.0 and archives any previous final folder.\n",
        "Step 04 alignment/orientation review decisions must be pass before assembly.\n", sep = "")
    quit(status = 0L)
  }
  config_env <- new.env()
  sys.source(find_wang_config(), envir = config_env)
  known <- c("validation-dir" = "validation_results_dir", "processing-dir" = "processing_input_dir",
             "overlay-dir" = "overlay_qa_dir", "release-dir" = "release_output_dir")
  seen <- character()
  for (arg in setdiff(raw_args, "--dry-run")) {
    require_release(grepl("^--[^=]+=.+$", arg), paste("Use --key=value:", arg))
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    require_release(key %in% names(known) && !key %in% seen, paste("Unknown or duplicate option:", key))
    seen <- c(seen, key)
    config_env$CONFIG[[known[[key]]]] <- sub("^--[^=]+=", "", arg)
  }
  release_results <- run_prepare_release(
    config_env$CONFIG, dry_run = "--dry-run" %in% commandArgs(TRUE)
  )
  if ("--dry-run" %in% commandArgs(TRUE) && !release_results$release_ready) {
    quit(status = 1L)
  }
}
