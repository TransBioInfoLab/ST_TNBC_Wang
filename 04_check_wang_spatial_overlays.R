#!/usr/bin/env Rscript
# Step 04: representative image/spot registration check BEFORE final assembly.
# Dependencies: EBImage (used to read the source Image objects); otherwise base R.
#
# From the ST_TNBC_Wang folder:
#   Rscript 04_check_wang_spatial_overlays.R
# Optional arguments (use --key=value):
#   --project-dir=/path/to/ST_TNBC_Wang
#   --processing-dir=/path/to/completed_processing_outputs
#   --out-dir=/path/to/a/NEW/qa_run_directory
#   --n=12
#   --arrays=CN1_C1,CN25_C2,CN43_C1,CN44_C1,CN44_C2
#
# Reads Step 02 outputs directly; no Step 05 release/manifest is required.
# Default input: outputs_wang_v1.1, or its existing Wang_TNBC_v1.1.0 subfolder.
# Default QA output: spatial_overlay_QA_v1.1.0/run_<UTC time>.
# Step 05 copies the chosen run into Wang_TNBC_v1.1.0/verification/spatial_overlays.
# --release-dir remains an alias for --processing-dir for previous commands.
# Existing output directories are refused, so earlier figures/review cannot be mixed.
# Source data and existing files are never modified. A custom QA destination
# inside a release must be under verification/spatial_overlays; Step 05 then reseals it.
#
# IMPORTANT: EBImage stores [x, y, channel], unlike the usual R [row, column].
# We transpose ONCE for raster display, retain the source pixel_x/pixel_y, and
# display y downwards. No extra rotation/reflection is applied to image or spots.
# The affine fit below DESCRIBES the author's grid-to-image mapping; it does not
# correct coordinates or establish biological registration from first principles.
# Full-resolution annotationRecoded masks are NOT overlaid: no documented mapping
# from their large raster frame to these per-array images has been established.
# The manual tumor panel uses existing per-SPOT annotations, not a raster boundary.
# This check does not approve pathology boundaries or validate cross-array distances.

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, ": missing ", paste(missing, collapse = ", "))
}

csv_out <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

full_path <- function(path) normalizePath(path, mustWork = TRUE)

payload_path <- function(root, relative) {
  if (length(relative) != 1L || is.na(relative) || !nzchar(relative) ||
      grepl("^(/|[A-Za-z]:)|(^|[/\\\\])\\.\\.([/\\\\]|$)", relative)) {
    stop("Invalid relative release path: ", relative)
  }
  path <- full_path(file.path(root, relative))
  if (!startsWith(path, paste0(root, "/"))) stop("Path escapes release: ", relative)
  path
}

checksum <- function(path) unname(tools::md5sum(path))

assert_sealed_input <- function(root, relative, manifest) {
  i <- match(relative, manifest$relative_path)
  if (is.na(i)) stop("Input is absent from the release checksum manifest: ", relative)
  path <- payload_path(root, relative)
  md5 <- checksum(path)
  if (is.na(md5) || !identical(md5, as.character(manifest$md5[i]))) {
    stop("Release input checksum mismatch: ", relative)
  }
  list(path = path, md5 = md5)
}

overlay_processing_root <- function(project_dir, requested = NULL) {
  if (!is.null(requested)) return(full_path(requested))
  config_file <- file.path(project_dir, "00_wang_config.R")
  cfg <- new.env(parent = baseenv())
  sys.source(config_file, envir = cfg)
  working <- cfg$CONFIG$working_output_dir
  if (!grepl("^/", working)) working <- file.path(project_dir, working)
  version <- sub("^v", "", cfg$CONFIG$release_version)
  if (grepl("^[0-9]+\\.[0-9]+$", version)) version <- paste0(version, ".0")
  candidates <- c(working, file.path(working, paste0("Wang_TNBC_v", version)))
  found <- candidates[file.exists(file.path(candidates, "logs", "processing_config.rds"))]
  if (!length(found)) stop("Completed Step 02 outputs not found; use --processing-dir=PATH")
  full_path(found[1])
}

spot_key <- function(x) paste(x$x, x$y, sep = "x")

image_shape <- function(image) {
  d <- dim(image)
  if (length(d) != 3L || d[3] != 3L) stop("Expected a three-channel EBImage H&E image")
  c(width = d[1], height = d[2])
}

image_raster <- function(image, lighten = 0) {
  # Transposition is essential: EBImage's FIRST dimension is horizontal x.
  a <- EBImage::imageData(image)
  if (any(!is.finite(a)) || min(a) < 0 || max(a) > 1) {
    stop("Expected finite H&E intensities between zero and one")
  }
  if (length(lighten) != 1L || !is.finite(lighten) || lighten < 0 || lighten > 1) {
    stop("lighten must be a scalar between zero and one")
  }
  # Display-only white blending improves contrast under the markers. Panel A
  # always retains the original image; no source pixels or coordinates change.
  if (lighten > 0) a <- a * (1 - lighten) + lighten
  grDevices::as.raster(aperm(a, c(2, 1, 3)))
}

open_png <- function(file, width, height, res) {
  # Use the platform default (Quartz on macOS). Forcing Cairo on macOS can
  # silently fail when the optional XQuartz libraries are not installed.
  before <- grDevices::dev.cur()
  grDevices::png(file, width = width, height = height, res = res)
  if (grDevices::dev.cur() == before) stop("PNG device did not open: ", file)
}

grid_mapping <- function(spots) {
  finite <- stats::complete.cases(spots[, c("new_x", "new_y", "pixel_x", "pixel_y")]) &
    is.finite(spots$new_x) & is.finite(spots$new_y) &
    is.finite(spots$pixel_x) & is.finite(spots$pixel_y)
  if (sum(finite) < 3L) stop("Too few finite source coordinates")
  X <- cbind(intercept = 1, new_x = spots$new_x[finite], new_y = spots$new_y[finite])
  if (qr(X)$rank < 3L) stop("Source spots do not span a two-dimensional grid")
  Y <- as.matrix(spots[finite, c("pixel_x", "pixel_y")])
  B <- qr.solve(X, Y)
  residual <- sqrt(rowSums((X %*% B - Y)^2))
  list(coefficients = B, determinant = det(B[2:3, , drop = FALSE]),
       angle = atan2(B[2, 2], B[2, 1]) * 180 / pi,
       rmse = sqrt(mean(residual^2)), max_residual = max(residual))
}

fraction_mean <- function(m, column, supported) {
  if (!column %in% names(m)) return(NA_real_)
  z <- m[[column]][supported & is.finite(m[[column]])]
  if (length(z)) mean(z) else NA_real_
}

inspect_array <- function(m, source, sid, expected_n) {
  s <- source$spots
  required <- c("x", "y", "new_x", "new_y", "pixel_x", "pixel_y")
  require_columns(s, required, paste(sid, "source image spots"))
  require_columns(m, c("spot_id", "subarray_id", required, "nCount_RNA", "qc_pass",
                      "manual_annotation_supported"), paste(sid, "released spots"))
  if (anyDuplicated(spot_key(s)) || anyDuplicated(m$spot_id) ||
      anyDuplicated(spot_key(m))) stop(sid, ": duplicate spot identifiers")
  if (anyNA(s[, c("x", "y")]) || anyNA(m[, c("x", "y")])) {
    stop(sid, ": missing lattice coordinates")
  }
  if (anyNA(m$subarray_id) || !all(m$subarray_id == sid) ||
      anyNA(m$spot_id) || !all(m$spot_id == paste0(sid, ".", spot_key(m)))) {
    stop(sid, ": inconsistent array/spot identifiers")
  }
  j <- match(spot_key(m), spot_key(s))
  exact_spot_set <- nrow(m) == nrow(s) && !anyNA(j) &&
    !anyNA(match(spot_key(s), spot_key(m))) && nrow(m) == expected_n
  # Missing/reordered source spots must never be compared by row position.
  delta <- if (!anyNA(j) && all(is.finite(as.matrix(m[, required]))) &&
               all(is.finite(as.matrix(s[j, required])))) {
    max(abs(as.matrix(m[, required]) - as.matrix(s[j, required])))
  } else Inf
  coordinate_match <- is.finite(delta) && delta <= 1e-6
  d <- image_shape(source$imgs)
  finite <- is.finite(m$pixel_x) & is.finite(m$pixel_y)
  inside <- finite & m$pixel_x >= 0.5 & m$pixel_x <= d["width"] + 0.5 &
    m$pixel_y >= 0.5 & m$pixel_y <= d["height"] + 0.5
  mapping <- grid_mapping(s)
  supported <- !is.na(m$manual_annotation_supported) & m$manual_annotation_supported
  valid_tumor <- if ("manual_fraction_Tumor" %in% names(m)) {
    sum(supported & is.finite(m$manual_fraction_Tumor))
  } else 0L
  technical <- exact_spot_set && coordinate_match && all(inside) &&
    is.finite(mapping$determinant) && abs(mapping$determinant) > 1e-8 &&
    mapping$max_residual <= 1
  data.frame(subarray_id = sid, n_spots = nrow(m), n_source_spots = nrow(s),
    image_width_px = unname(d["width"]), image_height_px = unname(d["height"]),
    source_spot_set_matches = exact_spot_set,
    exported_coordinates_match_source = coordinate_match,
    max_source_coordinate_difference = delta,
    n_nonfinite_pixel_coordinates = sum(!finite), n_spots_outside_image = sum(!inside & finite),
    grid_to_image_angle_deg = mapping$angle,
    reflected_relative_to_grid = mapping$determinant < 0,
    grid_to_image_determinant = mapping$determinant,
    affine_rmse_px = mapping$rmse, affine_max_residual_px = mapping$max_residual,
    n_supported_manual_spots = sum(supported), n_supported_tumor_fractions = valid_tumor,
    mean_supported_tumor_fraction = fraction_mean(m, "manual_fraction_Tumor", supported),
    mean_supported_necrosis_fraction = fraction_mean(m, "manual_fraction_Necrosis", supported),
    mean_whitespace_fraction = if ("manual_whitespace_fraction_all_pixels" %in% names(m) &&
                                  any(is.finite(m$manual_whitespace_fraction_all_pixels))) {
      mean(m$manual_whitespace_fraction_all_pixels, na.rm = TRUE)
    } else NA_real_,
    automatic_geometry_ok = technical, stringsAsFactors = FALSE)
}

choose_arrays <- function(inventory, n = 12L, requested = NULL) {
  x <- inventory[order(inventory$subarray_id), , drop = FALSE]
  selected <- character(); reasons <- character()
  add <- function(sid, reason) {
    if (!length(sid) || is.na(sid)) return(invisible(NULL))
    if (sid %in% selected) {
      i <- match(sid, selected); reasons[i] <<- paste(reasons[i], reason, sep = "; ")
    } else {
      selected <<- c(selected, sid); reasons <<- c(reasons, reason)
    }
  }
  if (!is.null(requested)) {
    if (anyDuplicated(requested) || any(!requested %in% x$subarray_id)) {
      stop("--arrays must contain unique, successfully processed array IDs")
    }
    for (sid in requested) add(sid, "User-selected representative array")
  } else {
    # Deterministic, outcome-blind selection; no clinical outcomes are loaded.
    # Include every reflected case in this cohort before filling other strata.
    for (sid in x$subarray_id[x$reflected_relative_to_grid]) add(sid, "Reflected source mapping")
    for (sid in x$subarray_id[!x$automatic_geometry_ok]) add(sid, "Automatic geometry flag")
    ordinary <- x[!x$reflected_relative_to_grid, , drop = FALSE]
    for (angle in c(0, 90, 180, -90)) {
      if (!nrow(ordinary)) break
      distance <- abs(((ordinary$grid_to_image_angle_deg - angle + 180) %% 360) - 180)
      # Prefer an annotated section within 15 degrees, so this stratum also
      # permits manual morphology review. Fall back to the nearest angle.
      candidates <- which(distance <= 15 & ordinary$n_supported_manual_spots > 0)
      i <- if (length(candidates)) {
        candidates[order(-ordinary$n_supported_manual_spots[candidates], distance[candidates],
                         ordinary$subarray_id[candidates])[1]]
      } else order(distance, ordinary$subarray_id)[1]
      add(ordinary$subarray_id[i], paste0("Source grid orientation near ", angle, " degrees"))
    }
    extreme <- function(column, decreasing, reason) {
      z <- x[is.finite(x[[column]]), , drop = FALSE]
      if (nrow(z)) add(z$subarray_id[order(if (decreasing) -z[[column]] else z[[column]],
                                          z$subarray_id)[1]], reason)
    }
    extreme("mean_supported_tumor_fraction", FALSE, "Low supported manual tumor fraction")
    extreme("mean_supported_tumor_fraction", TRUE, "High supported manual tumor fraction")
    extreme("mean_supported_necrosis_fraction", TRUE, "High manual necrosis fraction")
    extreme("mean_whitespace_fraction", TRUE, "High manual whitespace fraction")
    for (sid in x$subarray_id[order(x$n_spots, x$subarray_id)]) {
      if (length(selected) >= n) break
      add(sid, "Small spot footprint / additional coverage")
    }
    # The default cohort yields 12 cases. Explicit smaller --n is a preview only.
    selected <- head(selected, n); reasons <- head(reasons, n)
  }
  z <- x[match(selected, x$subarray_id), , drop = FALSE]
  z$selection_reason <- reasons
  z
}

plot_background <- function(raster, width, height, title, subtitle = "") {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, width + 0.5), ylim = c(height + 0.5, 0.5),
                        asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(raster, 0.5, height + 0.5, width + 0.5, 0.5, interpolate = FALSE)
  graphics::box(col = "grey65")
  graphics::axis(1, cex.axis = 0.6, col = "grey70")
  graphics::axis(2, cex.axis = 0.6, las = 1, col = "grey70")
  graphics::title(main = title, cex.main = 0.95, line = 1.2)
  graphics::mtext(subtitle, side = 3, line = 0.2, cex = 0.62)
}

legend_bottom <- function(...) {
  graphics::legend("bottom", inset = c(0, -0.17), xpd = NA, bty = "n", cex = 0.65, ...)
}

overlay_points <- function(m, ...) graphics::points(m$pixel_x, m$pixel_y, ...)

plot_spot_panel <- function(source, m, title, subtitle = "") {
  d <- image_shape(source$imgs)
  plot_background(image_raster(source$imgs, lighten = 0.35), d[1], d[2], title, subtitle)
  passed <- !is.na(m$qc_pass) & m$qc_pass
  overlay_points(m[passed, , drop = FALSE], pch = 21, cex = 1.00, lwd = 0.85,
                 col = "#111111", bg = "#00DEFF")
  overlay_points(m[!passed, , drop = FALSE], pch = 23, cex = 1.15, lwd = 1.10,
                 col = "#111111", bg = "#FF9C00")
  # Show increasing original lattice axes in the SOURCE image orientation.
  B <- grid_mapping(source$spots)$coefficients
  origin <- c(d[1] * 0.13, d[2] * 0.13)
  for (i in 2:3) {
    v <- B[i, ]; v <- v / sqrt(sum(v^2)) * min(d) * 0.07
    graphics::arrows(origin[1], origin[2], origin[1] + v[1], origin[2] + v[2],
                     length = 0.07, lwd = 1.6, col = "black")
    graphics::text(origin[1] + 1.3 * v[1], origin[2] + 1.3 * v[2],
                   labels = c("grid x+", "grid y+")[i - 1L], cex = 0.65, col = "black")
  }
  legend_bottom(legend = c("QC pass", "QC flagged"), pch = c(21, 23),
                col = "#111111", pt.bg = c("#00DEFF", "#FF9C00"), pt.cex = 1.1, horiz = TRUE)
}

plot_array <- function(source, m, info, file) {
  open_png(file, width = 2400, height = 2200, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 3.2, 3.5, 1.0), oma = c(1.3, 0, 3.2, 0))
  d <- image_shape(source$imgs); raster <- image_raster(source$imgs)
  overlay_raster <- image_raster(source$imgs, lighten = 0.35)
  plot_background(raster, d[1], d[2], "A. H&E", "Source image; x right, y down (pixels)")
  plot_spot_panel(source, m, "B. Released spot centers", "Lightened H&E; source coordinates without additional rotation/reflection")
  plot_background(overlay_raster, d[1], d[2], "C. Manual tumor annotation", "Per-spot Tumor fraction among biological annotation pixels; lightened H&E")
  good <- !is.na(m$manual_annotation_supported) & m$manual_annotation_supported
  if ("manual_fraction_Tumor" %in% names(m)) good <- good & is.finite(m$manual_fraction_Tumor) else good[] <- FALSE
  overlay_points(m[!good, , drop = FALSE], pch = 21, cex = 0.70, lwd = 0.70,
                 col = "#707070", bg = "#F5F5F5")
  palette <- grDevices::hcl.colors(101, "Blue-Red 3")
  if (any(good)) {
    z <- m$manual_fraction_Tumor[good]
    if (any(z < 0 | z > 1)) stop("Invalid manual tumor fractions")
    overlay_points(m[good, , drop = FALSE], pch = 21, cex = 0.95, lwd = 0.70,
                   col = "#202020", bg = palette[1L + round(z * 100)])
    legend_bottom(legend = c("0", "0.5", "1", "Unsupported / missing"),
                  pch = 21, col = c(rep("#202020", 3), "#707070"),
                  pt.bg = c(palette[c(1, 51, 101)], "#F5F5F5"), pt.cex = 1.1, horiz = TRUE)
  } else {
    graphics::mtext("No supported manual tumor annotations for this section", side = 1, line = 2.7, cex = 0.7)
  }
  plot_background(overlay_raster, d[1], d[2], "D. RNA library size", "UMI count per spot; logarithmic color scale; lightened H&E")
  z <- log10(1 + m$nCount_RNA); finite <- is.finite(z)
  palette <- grDevices::hcl.colors(101, "Viridis")
  if (any(finite)) {
    limits <- range(z[finite]); span <- diff(limits)
    index <- if (span > 0) 1L + round((z[finite] - limits[1]) / span * 100) else rep(51L, sum(finite))
    overlay_points(m[finite, , drop = FALSE], pch = 21, cex = 0.90, lwd = 0.65,
                   col = "#202020", bg = palette[index])
    ticks <- seq(limits[1], limits[2], length.out = 3)
    legend_bottom(legend = format(round(10^ticks - 1), big.mark = ",", trim = TRUE),
                  pch = 21, col = "#202020", pt.bg = palette[c(1, 51, 101)], pt.cex = 1.1, horiz = TRUE)
  }
  reflection <- if (info$reflected_relative_to_grid) "; reflected relative to grid" else ""
  graphics::mtext(sprintf("%s | %d spots | grid x+ angle %.1f degrees%s",
                          info$subarray_id, nrow(m), info$grid_to_image_angle_deg, reflection),
                  outer = TRUE, side = 3, line = 1.0, cex = 1.1, font = 2)
  graphics::mtext("Visual review pending. Spot annotations are not an approved tumor boundary. Distances remain within arrays.",
                  outer = TRUE, side = 1, line = 0.2, cex = 0.72)
}

run_wang_overlay_check <- function(project_dir, release_dir = NULL, out_dir = NULL,
                                    n = 12L, arrays = NULL, script_file = NULL,
                                    processing_dir = NULL) {
  if (!requireNamespace("EBImage", quietly = TRUE)) stop("Install EBImage before running this script")
  project_dir <- full_path(project_dir)
  if (!is.null(release_dir) && !is.null(processing_dir)) stop("Supply only one input-directory option")
  root <- overlay_processing_root(project_dir, if (!is.null(processing_dir)) processing_dir else release_dir)
  proc <- readRDS(file.path(root, "logs", "processing_config.rds"))
  version <- sub("^v", "", proc$release_version)
  if (grepl("^[0-9]+\\.[0-9]+$", version)) version <- paste0(version, ".0")
  if (is.null(out_dir)) out_dir <- file.path(project_dir, paste0("spatial_overlay_QA_v", version),
    paste0("run_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")))
  # Resolve the existing parent, including symlinks, before enforcing release isolation.
  if (dir.exists(out_dir) || file.exists(out_dir)) stop("Choose a NEW output directory: ", out_dir)
  ancestor <- dirname(path.expand(out_dir)); suffix <- basename(out_dir)
  while (!dir.exists(ancestor)) {
    suffix <- file.path(basename(ancestor), suffix); ancestor <- dirname(ancestor)
  }
  out_dir <- file.path(full_path(ancestor), suffix)
  if (out_dir == root || (startsWith(out_dir, paste0(root, "/")) &&
      !startsWith(out_dir, paste0(root, "/verification/spatial_overlays/")) &&
      out_dir != file.path(root, "verification", "spatial_overlays"))) {
    stop("QA inside the input directory must use verification/spatial_overlays (a new directory)")
  }
  if (length(n) != 1L || is.na(n) || n < 1 || n != as.integer(n)) stop("--n must be a positive integer")
  inputs <- list()
  remember <- function(role, path, md5 = checksum(path)) {
    inputs[[paste(role, path, sep = ":")]] <<- data.frame(role = role, path = full_path(path),
      file_size_bytes = file.info(path)$size, md5 = md5, stringsAsFactors = FALSE)
  }
  consumed <- function(relative, role) {
    path <- payload_path(root, relative)
    remember(role, path); path
  }
  remember("processing_config", file.path(root, "logs", "processing_config.rds"))
  summary <- read_csv(consumed("logs/processing_run_summary.csv", "processing_run_summary"))
  cohort <- read_csv(consumed("status/subarray_processing_status_latest.csv", "processing_status"))
  source_inputs <- read_csv(consumed("logs/processing_input_manifest.csv", "processing_input_manifest"))
  require_columns(summary, c("run_status", "full_cohort", "run_id", "n_successful", "n_spots"), "Step 02 run summary")
  require_columns(cohort, c("subarray_id", "count_id", "success", "n_spots",
                           "spot_table_rds_file", "run_id"), "Step 02 status")
  require_columns(source_inputs, c("input", "path", "md5"), "Step 02 inputs")
  if (nrow(summary) != 1L || !identical(summary$run_status, "complete") ||
      !isTRUE(summary$full_cohort) || !identical(summary$run_id, proc$run_id) ||
      anyNA(cohort$run_id) || !all(cohort$run_id == proc$run_id)) stop("A completed, coherent full-cohort Step 02 run is required")
  cohort <- cohort[!is.na(cohort$success) & cohort$success, , drop = FALSE]
  if (!nrow(cohort) || anyDuplicated(cohort$subarray_id)) stop("Empty or duplicated successful cohort")
  if (nrow(cohort) != summary$n_successful || sum(cohort$n_spots) != summary$n_spots) stop("Incoherent processing cohort totals")
  source_file <- function(role) {
    row <- source_inputs[source_inputs$input == role, , drop = FALSE]
    if (nrow(row) != 1L || !file.exists(row$path) || !identical(checksum(row$path), row$md5)) {
      stop("Source input differs from the frozen processing manifest: ", role)
    }
    row$path
  }
  ids_file <- source_file("ids")
  remember("source_array_lookup", ids_file)
  ids <- readRDS(ids_file)
  if (anyDuplicated(rownames(ids))) stop("Duplicated source lookup array IDs")
  require_columns(ids, c("id"), "Source array lookup")
  j <- match(cohort$subarray_id, rownames(ids))
  if (anyNA(j) || !all(as.character(ids$id[j]) == as.character(cohort$count_id))) {
    stop("Release array/specimen IDs do not match source image lookup")
  }
  if (!is.null(script_file)) remember("overlay_script", script_file)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE)
  writeLines("RUNNING: review is pending; do not treat partial output as completed QA.", file.path(out_dir, "RUN_STATUS.txt"))
  inventory <- list()
  parent_ids <- sort(unique(cohort$count_id))
  for (parent in parent_ids) {
    image_file <- source_file(paste0("images_", parent))
    remember("source_images", image_file)
    images <- readRDS(image_file)
    rows <- which(cohort$count_id == parent)
    for (i in rows) {
      sid <- cohort$subarray_id[i]
      if (!sid %in% names(images)) stop("Source image missing: ", sid)
      spot_file <- consumed(cohort$spot_table_rds_file[i], "processed_spot_metadata")
      m <- readRDS(spot_file)
      z <- inspect_array(m, images[[sid]], sid, cohort$n_spots[i])
      z$count_id <- parent
      z$spot_table_rds_file <- cohort$spot_table_rds_file[i]
      z$source_image_file <- full_path(image_file)
      inventory[[sid]] <- z
    }
    if (match(parent, parent_ids) %% 20L == 0L) message("Inspected image coordinates: ", length(inventory), " arrays")
  }
  inventory <- do.call(rbind, inventory); rownames(inventory) <- NULL
  selected <- choose_arrays(inventory, n, arrays)
  selected$figure_file <- paste0("figures/", selected$subarray_id, "_overlay.png")
  csv_out(inventory, file.path(out_dir, "automatic_coordinate_checks_all_arrays.csv"))
  csv_out(selected, file.path(out_dir, "representative_selection.csv"))
  for (i in seq_len(nrow(selected))) {
    info <- selected[i, , drop = FALSE]
    images <- readRDS(info$source_image_file)
    m <- readRDS(payload_path(root, info$spot_table_rds_file))
    plot_array(images[[info$subarray_id]], m, info, file.path(out_dir, info$figure_file))
    message("Rendered ", info$subarray_id, " (", i, "/", nrow(selected), ")")
  }
  # A compact overview; detailed review should use the individual four-panel figures.
  ncol <- min(4L, nrow(selected)); nrow_plot <- ceiling(nrow(selected) / ncol)
  contact_file <- file.path(out_dir, "overlay_contact_sheet.png")
  open_png(contact_file, width = ncol * 1000, height = nrow_plot * 1000, res = 130)
  tryCatch({
    graphics::par(mfrow = c(nrow_plot, ncol), mar = c(3.6, 2.6, 3, 0.5), oma = c(0, 0, 2, 0))
    for (i in seq_len(nrow(selected))) {
      info <- selected[i, , drop = FALSE]; images <- readRDS(info$source_image_file)
      m <- readRDS(payload_path(root, info$spot_table_rds_file))
      plot_spot_panel(images[[info$subarray_id]], m, info$subarray_id,
        sprintf("%.1f degrees%s", info$grid_to_image_angle_deg,
          if (info$reflected_relative_to_grid) "; reflected mapping" else ""))
    }
    graphics::mtext("Wang representative image/spot overlays - VISUAL REVIEW PENDING", outer = TRUE,
                    side = 3, line = 0.3, cex = 1.0)
  }, finally = grDevices::dev.off())
  figure_paths <- c(file.path(out_dir, selected$figure_file), contact_file)
  if (any(!file.exists(figure_paths)) || anyNA(file.info(figure_paths)$size) ||
      any(file.info(figure_paths)$size <= 0)) stop("Missing or empty PNG output; this QA run is incomplete")
  review <- selected[, c("subarray_id", "count_id", "selection_reason", "figure_file",
                         "reflected_relative_to_grid", "automatic_geometry_ok")]
  review$image_alignment_review <- "pending"
  review$rotation_flip_review <- "pending"
  review$manual_annotation_review <- ifelse(selected$n_supported_tumor_fractions > 0, "pending", "not_available")
  csv_out(review, file.path(out_dir, "overlay_review.csv"))
  inputs <- do.call(rbind, inputs); rownames(inputs) <- NULL
  # Ensure nothing used for scanning, selection, or rendering changed during this run.
  unchanged <- unname(tools::md5sum(inputs$path)) == inputs$md5
  if (anyNA(unchanged) || !all(unchanged)) stop("Input changed during QA; do not use this run")
  csv_out(inputs, file.path(out_dir, "overlay_input_manifest.csv"))
  if (!is.null(script_file)) file.copy(script_file, file.path(out_dir, "overlay_script.R"), overwrite = FALSE)
  capture.output(utils::sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
  run_summary <- data.frame(release_version = version,
    processing_run_id = proc$run_id, processing_directory = root,
    input_mode = "completed_step02_outputs",
    run_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    n_arrays_checked = nrow(inventory), n_arrays_with_automatic_flags = sum(!inventory$automatic_geometry_ok),
    n_representative_arrays = nrow(selected),
    n_reflected_arrays_in_cohort = sum(inventory$reflected_relative_to_grid),
    n_reflected_arrays_selected = sum(selected$reflected_relative_to_grid),
    all_reflected_cases_selected = all(inventory$subarray_id[inventory$reflected_relative_to_grid] %in% selected$subarray_id),
    n_selected_with_supported_manual_tumor = sum(selected$n_supported_tumor_fractions > 0),
    visual_review_status = "pending", approved_pathology_boundary = FALSE,
    cross_array_registration_validated = FALSE, stringsAsFactors = FALSE)
  csv_out(run_summary, file.path(out_dir, "qa_run_summary.csv"))
  writeLines(c(
    "Step 04: Wang representative image/spot overlay QA", "",
    "Read completed Step 02 outputs directly; final Step 05 assembly is not needed to run this check.",
    "Panels B-D use larger opaque markers and H&E lightened by 35%; panel A retains the source image.",
    "Marker size is for visibility, not a calibrated capture radius. Source data/coordinates are unchanged.",
    "Review figures and fill overlay_review.csv: image_alignment_review and rotation_flip_review",
    "Grey/white manual-annotation points are unsupported/missing, not tumor-negative. Only the Tumor category",
    "is shown, restricted to supported spots. This does not approve full-resolution masks or invasive boundaries.",
    "The coordinate audit covers all successful arrays; representative figures cover the selected arrays only.",
    "Source mapping rotations/reflections are retained; distances remain within arrays.", "",
    "Next: Rscript 05_prepare_wang_release.R --overlay-dir=/absolute/path/to/this/run --dry-run",
    "Then run the same command without --dry-run after the visual decisions have been recorded.",
    "Step 05 copies only this chosen run into Wang_TNBC_v1.1.0/verification/spatial_overlays, updates release",
    "documentation and summaries, and verifies a fresh full inventory/checksum manifest. Existing releases are",
    "rebuilt in staging and retained as a backup before the original folder name is restored.",
    "qa_output_manifest.csv protects generated output except the editable overlay_review.csv and the manifest",
    "itself. Step 05 separately captures and validates the final review CSV. Do not edit protected figures/tables;",
    "generate a new QA run instead. Reruns use fresh directories to avoid mixing old and new output."
  ), file.path(out_dir, "README_overlay_check.txt"))
  writeLines("COMPLETE: automatic checks/figures generated; human visual review is PENDING.", file.path(out_dir, "RUN_STATUS.txt"))
  payload <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
  payload <- sort(setdiff(payload, c("overlay_review.csv", "qa_output_manifest.csv")))
  csv_out(data.frame(relative_path = payload, file_size_bytes = file.info(file.path(out_dir, payload))$size,
                     md5 = unname(tools::md5sum(file.path(out_dir, payload)))),
          file.path(out_dir, "qa_output_manifest.csv"))
  message("QA output: ", out_dir)
  message("Automatic geometry flags: ", sum(!inventory$automatic_geometry_ok), "/", nrow(inventory),
          ". Representative figures: ", nrow(selected), ". Visual review remains PENDING.")
  invisible(list(output_dir = out_dir, inventory = inventory, selected = selected, summary = run_summary))
}

if (sys.nframe() == 0L) {
  raw_args <- commandArgs(trailingOnly = TRUE)
  if (any(raw_args %in% c("--help", "-h"))) {
    cat("Usage: Rscript 04_check_wang_spatial_overlays.R [--project-dir=PATH] [--processing-dir=PATH]\n",
        "       [--out-dir=NEW_PATH] [--n=12] [--arrays=CN1_C1,CN25_C2,...]\n",
        "Reads Step 02 output before release assembly. Step 05 incorporates the selected QA run.\n", sep = "")
  } else {
    args <- list()
    for (a in raw_args) {
      if (!grepl("^--[^=]+=.+$", a)) stop("Use --key=value: ", a)
      key <- sub("^--([^=]+)=.*$", "\\1", a); value <- sub("^--[^=]+=", "", a)
      if (!key %in% c("project-dir", "processing-dir", "release-dir", "out-dir", "n", "arrays") || key %in% names(args)) {
        stop("Unknown or duplicate option: ", key)
      }
      args[[key]] <- value
    }
    script_arg <- grep("^--file=", commandArgs(), value = TRUE)
    script_file <- if (length(script_arg)) full_path(sub("^--file=", "", script_arg[1])) else NULL
    project <- if (!is.null(args[["project-dir"]])) args[["project-dir"]] else {
      if (!is.null(script_file)) dirname(script_file) else getwd()
    }
    run_wang_overlay_check(project, release_dir = args[["release-dir"]], out_dir = args[["out-dir"]],
      n = if (is.null(args$n)) 12L else as.numeric(args$n),
      arrays = if (is.null(args$arrays)) NULL else strsplit(args$arrays, ",", fixed = TRUE)[[1]],
      script_file = script_file, processing_dir = args[["processing-dir"]])
  }
}
