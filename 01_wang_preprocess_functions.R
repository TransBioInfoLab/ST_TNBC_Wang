################################################################################
# TNBC spatial transcriptomics preprocessing in Wang's paper

# set library
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})



# ---- helper functions ----

# create key spot
make_wang_spot_key <- function(slide, x, y) {
  paste0(slide, ".", x, "x", y)
}

# load parent object
load_wang_parent <- function(
    count_id,
    ids,
    clinical,
    counts_dir,
    images_dir,
    classification,
    annotation_dir = NULL
) {
  
  count_id <- as.character(count_id)
  count_file <- file.path(counts_dir, paste0("TNBC", count_id, ".RDS"))
  image_file <- file.path(images_dir, paste0("TNBC", count_id, ".RDS"))
  
  if (!file.exists(count_file)) {
    stop("Count file not found: ", count_file)
  }
  
  if (!file.exists(image_file)) {
    stop("Image file not found: ", image_file)
  }
  
  count_obj <- readRDS(count_file)
  image_obj <- readRDS(image_file)
  
  manifest_rows <- ids[as.character(ids$id) == count_id, ,drop = FALSE]
  
  if (nrow(manifest_rows) == 0) {
    stop("No ids rows found for count ID: ", count_id)
  }
  
  clinical_index <- match(
    count_id,
    as.character(clinical$ST_TNBC_ID)
  )
  
  clinical_row <- if (!is.na(clinical_index)) {
    clinical[clinical_index, , drop = FALSE]
  } else {
    NULL
  }
  
  class_obj <- if (count_id %in% names(classification)) {
    classification[[count_id]]
  } else {
    NULL
  }
  
  annotation_obj <- NULL
  annotation_file <- file.path(annotation_dir, paste0("TNBC", count_id, ".RDS"))
  
  if (!is.null(annotation_dir) &&file.exists(annotation_file)) {
    annotation_obj <- readRDS(annotation_file)
  }
  
  # Construct expression spot keys
  count_obj$spots$spot_key <- paste0(
    count_obj$spots$slide, ".", count_obj$spots$x, "x", count_obj$spots$y
  )
  
  matrix_spot_keys <- paste0(
    count_obj$spots$slide,
    ".",
    rownames(count_obj$cnts)
  )
  
  if (!identical(matrix_spot_keys, count_obj$spots$spot_key)) {
    stop(
      "Count matrix rows and count_obj$spots do not align ",
      "for count ID: ", count_id
    )
  }
  
  # Attach classification identifiers
  if (!is.null(class_obj)) {
    
    if (
      nrow(class_obj$pr) != nrow(class_obj$spot) ||
      nrow(class_obj$pr) != nrow(class_obj$xy)
    ) {
      stop(
        "Classification dimensions are inconsistent for ",
        "count ID: ", count_id
      )
    }
    
    class_keys <- rownames(class_obj$spot)
    rownames(class_obj$pr) <- class_keys
    rownames(class_obj$xy) <- class_keys
  }
  
  list(
    count_id = count_id,
    count_obj = count_obj,
    image_obj = image_obj,
    classification_obj = class_obj,
    annotation_obj = annotation_obj,
    manifest_rows = manifest_rows,
    clinical_row = clinical_row
  )
}


# Convert patient-level clinical values into ordinary scalar metadata columns.
# survival::Surv objects are split into separate time/status columns so that
# spot metadata remains a regular rectangular data frame.
add_clinical_metadata_safe <- function(obj, clinical_row) {
  if (is.null(clinical_row)) return(obj)

  for (nm in names(clinical_row)) {
    value <- clinical_row[[nm]]

    if (inherits(value, "Surv")) {
      value_matrix <- unclass(value)
      value_names <- colnames(value_matrix)

      if (is.null(value_names)) {
        value_names <- paste0("value", seq_len(ncol(value_matrix)))
      }

      for (j in seq_len(ncol(value_matrix))) {
        out_name <- paste0(
          "clinical_", nm, "_",
          make.names(value_names[j], unique = TRUE)
        )

        obj[[out_name]] <- rep(as.numeric(value_matrix[1, j]), ncol(obj))
      }

      next
    }

    flat <- unlist(value, recursive = TRUE, use.names = FALSE)

    if (length(flat) == 0) {
      flat <- NA_character_
    } else if (length(flat) > 1) {
      flat <- paste(as.character(flat), collapse = ";")
    }

    obj[[paste0("clinical_", nm)]] <- rep(flat[[1]], ncol(obj))
  }

  obj
}

# load subarray
extract_wang_subarray <- function(parent, subarray_id){
  
  count_obj <- parent$count_obj
  image_obj <- parent$image_obj
  class_obj <- parent$classification_obj
  annotation_obj <- parent$annotation_obj
  manifest_rows <- parent$manifest_rows
  clinical_row <- parent$clinical_row
  count_id <- parent$count_id
  
  if (!subarray_id %in% rownames(manifest_rows)) {
    stop(
      "Subarray ", subarray_id,
      " does not belong to count ID ", count_id
    )
  }
  
  if (!subarray_id %in% names(image_obj)) {
    stop(
      "No image found for subarray: ",
      subarray_id
    )
  }
  
  id_row <- manifest_rows[subarray_id, ,drop = FALSE]
  keep <- count_obj$spots$slide == subarray_id
  
  if (!any(keep)) {
    stop(
      "No spots found for subarray: ",
      subarray_id
    )
  }
  
  spot_meta <- count_obj$spots[keep, ,drop = FALSE]
  counts_spot_gene <- count_obj$cnts[keep, ,drop = FALSE]
  
  rownames(counts_spot_gene) <- spot_meta$spot_key
  rownames(spot_meta) <- spot_meta$spot_key
  
  # Seurat expects genes × spots
  counts_gene_spot <- t(counts_spot_gene)
  
  # Convert to sparse format.
  counts_gene_spot <- Matrix::Matrix(counts_gene_spot, sparse = TRUE)
  detected_genes <- Matrix::rowSums(counts_gene_spot > 0) > 0
  counts_gene_spot <- counts_gene_spot[detected_genes, ,drop = FALSE]
  
  obj <- Seurat::CreateSeuratObject(
    counts = counts_gene_spot,
    project = subarray_id,
    min.cells = 1,
    min.features = 0,
    meta.data = spot_meta
  )
  
  obj$study <- "Wang_2024"
  obj$platform <- "legacy_ST"
  obj$count_id <- count_id
  obj$subarray_id <- subarray_id
  obj$parent_slide <- as.character(id_row$slide)
  obj$subslide <- as.character(id_row$subslide)
  obj$sample_label <- as.character(id_row$ID)
  obj$has_manual_annotation <-
    as.logical(id_row$hasAnnot)
  
  # Add clinical metadata in scalar-safe form.
  obj <- add_clinical_metadata_safe(
    obj = obj,
    clinical_row = clinical_row
  )

  # Add classification scores
  if (!is.null(class_obj)) {
    
    class_index <- match(colnames(obj), rownames(class_obj$pr))
    
    if (anyNA(class_index)) {
      stop(
        sum(is.na(class_index)),
        " spots did not match classification scores ",
        "for subarray ", subarray_id
      )
    }
    
    score_matrix <- class_obj$pr[class_index, ,drop = FALSE]
    original_names <- colnames(score_matrix)
    colnames(score_matrix) <- paste0(
      "hist_score_",
      make.names(original_names, unique = TRUE)
    )
    
    score_df <- as.data.frame(score_matrix)
    rownames(score_df) <- colnames(obj)
    obj <- AddMetaData(obj,metadata = score_df)
    
    obj$hist_highest_score_class <- original_names[max.col(score_matrix, ties.method = "first")]
    obj$hist_highest_score_value <- apply(score_matrix, 1, max, na.rm = TRUE)
    obj$hist_score_sum <- rowSums(score_matrix, na.rm = TRUE)
  }
  
  
  # Add manual histopathology annotations
  if (isTRUE(as.logical(id_row$hasAnnot))) {
    
    if (is.null(annotation_obj)) {
      warning(
        "ids indicates manual annotation is available for ",
        subarray_id,
        ", but no annotation object was loaded."
      )
      
    } else {
      
      if (!all(c("annots", "spots") %in% names(annotation_obj))) {
        stop(
          "Annotation object for count ID ",
          count_id,
          " must contain both `annots` and `spots`."
        )
      }
      
      if (nrow(annotation_obj$annots) != nrow(annotation_obj$spots)) {
        stop(
          "Annotation matrix and annotation spot metadata ",
          "have different row counts for count ID: ",
          count_id
        )
      }
      
      # Construct unique annotation spot keys
      annotation_spot_keys <- paste0(subarray_id, ".",
        annotation_obj$spots$x, "x", annotation_obj$spots$y)
      
      annotation_matrix <- annotation_obj$annots
      rownames(annotation_matrix) <- annotation_spot_keys
      
      # Match manual annotations to the current Seurat object
      annotation_index <- match(colnames(obj), annotation_spot_keys)
      
      if (anyNA(annotation_index)) {
        
        missing_spots <- colnames(obj)[is.na(annotation_index)]
        
        stop(length(missing_spots), " spots in ", subarray_id,
          " could not be matched to manual annotations. ",
          "Examples: ",paste(head(missing_spots, 5), collapse = ", "))
      }
      
      annotation_counts <- annotation_matrix[annotation_index, ,drop = FALSE]
      rownames(annotation_counts) <- colnames(obj)
      
      # define excluded annotation categories
      excluded_categories <- c("Nothing", "Hole (whitespace)","Artefacts")
      biological_categories <- setdiff(colnames(annotation_counts), excluded_categories)
      
      total_all_pixels <- rowSums(annotation_counts,na.rm = TRUE)
      total_biological_pixels <- rowSums(
        annotation_counts[, biological_categories, drop = FALSE],
        na.rm = TRUE
      )
      
      # calculate biological fractions
      biological_counts <- annotation_counts[, biological_categories, drop = FALSE]
      biological_fraction <- sweep(biological_counts, MARGIN = 1, STATS = pmax(total_biological_pixels, 1), FUN = "/")
      
      # Keep the original category names before sanitizing them
      original_annotation_names <- colnames(annotation_counts)
      safe_annotation_names <- make.names(original_annotation_names, unique = TRUE)
      
      # Add raw annotation pixel counts
      annotation_counts_df <- as.data.frame(annotation_counts)
      colnames(annotation_counts_df) <- paste0("manual_pixels_", safe_annotation_names)
      rownames(annotation_counts_df) <- colnames(obj)
      
      obj <- Seurat::AddMetaData(object = obj, metadata = annotation_counts_df)
      
      # Add within-spot annotation fractions
      biological_safe_names <- make.names(biological_categories, unique = TRUE)
      annotation_fraction_df <- as.data.frame(biological_fraction)
      colnames(annotation_fraction_df) <- paste0("manual_fraction_", biological_safe_names)
      rownames(annotation_fraction_df) <- colnames(obj)
      
      obj <- Seurat::AddMetaData(object = obj, metadata = annotation_fraction_df)
      
      # calculte QC fractions
      artifact_fraction <- annotation_counts[, "Artefacts"] / pmax(total_all_pixels, 1)
      whitespace_fraction <- annotation_counts[, "Hole (whitespace)"] / pmax(total_all_pixels, 1)
      unannotated_fraction <- annotation_counts[, "Nothing"] / pmax(total_all_pixels, 1)
      
      obj$manual_artifact_fraction_all_pixels <- artifact_fraction
      obj$manual_whitespace_fraction_all_pixels <- whitespace_fraction
      obj$manual_unannotated_fraction_all_pixels <- unannotated_fraction
      
      # Dominant manually annotated category
      manual_dominant_class <- rep(NA_character_, nrow(biological_counts))
      annotated <- total_biological_pixels > 0
      
      manual_dominant_class[annotated] <-
        colnames(biological_counts)[
          max.col(biological_counts[annotated, ,drop = FALSE], ties.method = "first")]
      
      obj$manual_dominant_class <- manual_dominant_class
      obj$manual_all_pixels <- total_all_pixels
      obj$manual_annotated_pixels <- total_biological_pixels
      obj$manual_has_annotated_pixels <- total_biological_pixels > 0
    }
    
  } else {
    
    # This subarray was not manually annotated.
    # Add explicit status variables for consistency.
    obj$manual_artifact_fraction_all_pixels <- NA_real_
    obj$manual_whitespace_fraction_all_pixels <- NA_real_
    obj$manual_unannotated_fraction_all_pixels <- NA_real_
    obj$manual_dominant_class <- NA_character_
    obj$manual_all_pixels <- NA_real_
    obj$manual_annotated_pixels <- NA_real_
    obj$manual_has_annotated_pixels <- FALSE
  }
  
  list(
    seurat = obj,
    image = image_obj[[subarray_id]]$imgs,
    image_spots = image_obj[[subarray_id]]$spots,
    sample_info = id_row
  )
}

# Wang's count matrics are already filtered, so just record the qc
perform_wang_qc <- function(
    obj,
    min_umi = 500,
    max_artifact_fraction = 0.01
) {
  
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj,pattern = "^MT-")
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj,pattern = "^RP[SL]")
  obj[["percent.mtrnr"]] <- Seurat::PercentageFeatureSet(obj,pattern = "^MTRNR")
  
  obj$qc_pass_umi <- obj$nCount_RNA >= min_umi
  
  if ("manual_artifact_fraction_all_pixels" %in% colnames(obj@meta.data)) {
    obj$qc_pass_artifact <-
      is.na(obj$manual_artifact_fraction_all_pixels) |
      obj$manual_artifact_fraction_all_pixels <= max_artifact_fraction
  } else {
    obj$qc_pass_artifact <- TRUE
  }
  
  obj$qc_pass <- obj$qc_pass_umi & obj$qc_pass_artifact
  
  obj
}

# convert coordinates to microns
add_wang_micron_coordinates <- function(obj,
                                        spot_pitch_um = 150) {
  
  # instead deriving a different scale for every subarray, used a fixed constant
  um_per_unit <- spot_pitch_um / sqrt(2)
  
  obj$x_um <- (
    obj$x - min(obj$x, na.rm = TRUE)
  ) * um_per_unit
  
  obj$y_um <- (
    obj$y - min(obj$y, na.rm = TRUE)
  ) * um_per_unit
  
  obj$coordinate_source <- "original_array_xy"
  obj$microns_per_coordinate_unit <- um_per_unit
  obj$assumed_spot_pitch_um <- spot_pitch_um
  
  obj
}

# normalize each subarray
normalize_wang_subarray <- function(obj) {
  
  obj <- Seurat::NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  obj
}

# compute module scores
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

get_assay_matrix <- function(obj, assay = NULL, layer = "data") {
  assay <- assay %||% DefaultAssay(obj)
  # Compatible with Seurat v4/v5. LayerData exists in SeuratObject >=5.
  if ("LayerData" %in% getNamespaceExports("SeuratObject")) {
    mat <- tryCatch(
      SeuratObject::LayerData(obj, assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (!is.null(mat)) return(mat)
  }
  # Fallback for older Seurat APIs.
  suppressWarnings({
    mat <- Seurat::GetAssayData(obj, assay = assay, slot = layer)
  })
  mat
}

robust_z <- function(x) {
  x <- as.numeric(x)
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(madv) || madv <= 0) {
    sdv <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sdv) || sdv <= 0) return(rep(0, length(x)))
    return((x - mean(x, na.rm = TRUE)) / sdv)
  }
  (x - med) / madv
}

nonnegative_scale <- function(x) {
  # Convert any abundance-like score to a nonnegative field.
  # This is a robust z-score truncated at zero; it is not a cell proportion.
  y <- robust_z(x)
  y[!is.finite(y)] <- 0
  pmax(y, 0)
}

calculate_module_coverage <- function(obj, gene_sets, subarray_id, count_id) {
  data.table::rbindlist(lapply(names(gene_sets), function(nm) {
    requested <- unique(gene_sets[[nm]])
    detected <- intersect(requested, rownames(obj))
    data.frame(
      count_id = as.character(count_id),
      subarray_id = as.character(subarray_id),
      module = nm,
      n_requested = length(requested),
      n_detected = length(detected),
      fraction_detected = if (length(requested) > 0) length(detected) / length(requested) else NA_real_,
      missing_genes = paste(setdiff(requested, detected), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)
}

compute_module_scores <- function(obj, gene_sets, assay = NULL, layer = "data") {
  assay <- assay %||% DefaultAssay(obj)
  expr <- get_assay_matrix(obj, assay = assay, layer = layer)
  all_genes <- rownames(expr)
  
  meta_add <- list()
  
  for (nm in names(gene_sets)) {
    genes <- intersect(gene_sets[[nm]], all_genes)
    out_col <- paste0("score_", nm)
    
    if (length(genes) == 0) {
      warning("No genes found for module: ", nm)
      meta_add[[out_col]] <- rep(NA_real_, ncol(expr))
    } else if (length(genes) == 1) {
      meta_add[[out_col]] <- as.numeric(expr[genes, ])
    } else {
      meta_add[[out_col]] <- Matrix::colMeans(expr[genes, , drop = FALSE])
    }
  }
  
  meta_df <- as.data.frame(meta_add)
  rownames(meta_df) <- colnames(expr)
  
  # Add robust z-score version for every module score.
  for (nm in colnames(meta_df)) {
    meta_df[[paste0(nm, "_z")]] <- robust_z(meta_df[[nm]])
    meta_df[[paste0(nm, "_pos")]] <- nonnegative_scale(meta_df[[nm]])
  }
  
  Seurat::AddMetaData(obj, metadata = meta_df)
}

# Ensure a spot-level table contains exactly one scalar value per spot.
validate_scalar_spot_table <- function(df, table_name = "spot table") {
  df <- as.data.frame(df, check.names = FALSE)
  nr <- nrow(df)

  bad <- vapply(
    df,
    function(x) {
      !is.atomic(x) ||
        !is.null(dim(x)) ||
        length(x) != nr
    },
    logical(1)
  )

  if (any(bad)) {
    details <- vapply(
      names(df)[bad],
      function(nm) {
        x <- df[[nm]]
        paste0(
          nm,
          " [class=", paste(class(x), collapse = "/"),
          ", length=", length(x),
          if (!is.null(dim(x))) {
            paste0(", dim=", paste(dim(x), collapse = "x"))
          } else {
            ""
          },
          "]"
        )
      },
      character(1)
    )

    stop(
      table_name,
      " contains non-scalar columns: ",
      paste(details, collapse = "; ")
    )
  }

  invisible(TRUE)
}

