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
extract_wang_subarray <- function(parent, subarray_id, 
                                  min_annotation_pixels = 1000,
                                  relapse_to_primary = c("30" = "53", "66" = "58")){
  
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
  
  # Seurat replaces underscores and pipes in feature names.
  # Record the conversion while preserving expression values.
  source_genes <- rownames(counts_gene_spot)
  if (is.null(source_genes) || anyNA(source_genes) ||
      any(!nzchar(source_genes)) || anyDuplicated(source_genes)) {
    stop("Source feature names must be present and unique.")
  }
  
  assay_genes <- gsub("[_|]", "-", source_genes)
  if (anyDuplicated(assay_genes)) {
    stop("Feature names collide after Seurat name conversion.")
  }
  
  rownames(counts_gene_spot) <- assay_genes
  obj <- Seurat::CreateSeuratObject(
    counts = counts_gene_spot,
    project = subarray_id,
    min.cells = 0,
    min.features = 0,
    meta.data = spot_meta
  )
  
  obj@misc$wang_feature_map <- data.frame(
    source_gene = source_genes,
    assay_gene = assay_genes,
    stringsAsFactors = FALSE
  )
  
  obj$study <- "Wang_2024"
  obj$platform <- "legacy_ST"
  obj$count_id <- count_id
  obj$subarray_id <- subarray_id
  
  is_relapse <- count_id %in% names(relapse_to_primary)
  obj$subject_id <- if (is_relapse) {
    unname(relapse_to_primary[count_id])
  } else {count_id}
  
  obj$specimen_type <- if (is_relapse) "relapse" else "primary"
  obj$baseline_primary_specimen <- !is_relapse
  
  obj$parent_slide <- as.character(id_row$slide)
  obj$subslide <- as.character(id_row$subslide)
  obj$sample_label <- as.character(id_row$ID)
  obj$has_manual_annotation <- as.logical(id_row$hasAnnot)
  
  # Add clinical metadata in scalar-safe form.
  obj <- add_clinical_metadata_safe(obj = obj, clinical_row = clinical_row)

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
    # These labels describe the current classifAll.RDS input.
    obj$hist_prediction_source <- "classifAll.RDS"
    obj$hist_prediction_type <- "full_data_expression_prediction"
    obj$hist_prediction_modality <- "RNA"
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
      biological_counts <- annotation_counts[, biological_categories, drop = FALSE]
      
      total_all_pixels <- rowSums(annotation_counts,na.rm = TRUE)
      total_biological_pixels <- rowSums(biological_counts)
      
      # calculate biological fractions
      bio_denominator <- total_biological_pixels
      bio_denominator[!is.finite(bio_denominator) | bio_denominator <= 0] <- NA_real_
      
      biological_fraction <- sweep(biological_counts, 1, bio_denominator, "/")
      
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
      all_denominator <- total_all_pixels
      all_denominator[!is.finite(all_denominator) | all_denominator <= 0] <- NA_real_
      
      artifact_fraction <- annotation_counts[, "Artefacts"] / all_denominator
      whitespace_fraction <- annotation_counts[, "Hole (whitespace)"] / all_denominator
      unannotated_fraction <- annotation_counts[, "Nothing"] / all_denominator
      
      obj$manual_artifact_fraction_all_pixels <- artifact_fraction
      obj$manual_whitespace_fraction_all_pixels <- whitespace_fraction
      obj$manual_unannotated_fraction_all_pixels <- unannotated_fraction
      
      # Dominant manually annotated category
      manual_dominant_class <- rep(NA_character_, nrow(biological_counts))
      annotated <- is.finite(total_biological_pixels) & total_biological_pixels > 0
      
      manual_dominant_class[annotated] <-
        colnames(biological_counts)[
          max.col(biological_counts[annotated, ,drop = FALSE], ties.method = "first")]
      
      obj$manual_dominant_class <- manual_dominant_class
      obj$manual_all_pixels <- total_all_pixels
      obj$manual_annotated_pixels <- total_biological_pixels
      
      obj$manual_has_annotated_pixels <- annotated
      obj$manual_annotation_supported <- is.finite(total_biological_pixels) &
        total_biological_pixels > min_annotation_pixels
      obj$manual_annotation_min_pixels <- min_annotation_pixels
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
    obj$manual_annotation_supported <- FALSE
    obj$manual_annotation_min_pixels <- min_annotation_pixels
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
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  
  if (!any(ok)) return(out)
  
  values <- x[ok]
  med <- stats::median(values)
  spread <- stats::mad(values, center = med, constant = 1.4826)
  
  if (is.finite(spread) && spread > 0) {
    out[ok] <- (values - med) / spread
  } else {
    spread <- stats::sd(values)
    
    out[ok] <- if (is.finite(spread) && spread > 0) {
      (values - mean(values)) / spread
    } else {
      rep(0, length(values))
    }
  }
  
  out
}

nonnegative_scale <- function(x) {
  pmax(robust_z(x), 0)
}

wang_gene_alias_groups <- function() {
  list(
    NIBAN3 = c("NIBAN3", "FAM129C"),
    RIPOR2 = c("RIPOR2", "FAM65B"),
    ACKR3 = c("ACKR3", "CXCR7")
  )
}

resolve_wang_genes <- function(requested, feature_map,
    alias_groups = wang_gene_alias_groups()) {
  requested <- unique(as.character(requested))
  
  if (!length(requested) || anyNA(requested) ||
      any(!nzchar(requested))) {
    stop("Each module must contain nonempty gene identifiers.")
  }
  
  if (!all(c("source_gene", "assay_gene") %in% names(feature_map)) ||
      anyNA(feature_map) ||
      anyDuplicated(feature_map$source_gene) ||
      anyDuplicated(feature_map$assay_gene)) {
    stop("Invalid source-to-assay feature map.")
  }
  
  source_gene <- vapply(requested, function(g) {
    groups <- alias_groups[vapply(
      alias_groups, function(a) g %in% a, logical(1)
    )]
    
    candidates <- unique(c(
      g, unlist(groups, use.names = FALSE)
    ))
    
    hits <- intersect(candidates, feature_map$source_gene)
    
    if (length(hits) > 1L) {
      stop(
        "Ambiguous mapping for ", g, ": ",
        paste(hits, collapse = ", ")
      )
    }
    
    if (length(hits)) hits[[1L]] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  
  assay_gene <- feature_map$assay_gene[
    match(source_gene, feature_map$source_gene)
  ]
  
  if (anyDuplicated(assay_gene[!is.na(assay_gene)])) {
    stop(
      "A module maps multiple requested genes ",
      "to the same assay feature."
    )
  }
  
  data.frame(
    requested_gene = requested,
    source_gene = source_gene,
    assay_gene = assay_gene,
    alias_used = !is.na(source_gene) &
      requested != source_gene,
    stringsAsFactors = FALSE
  )
}

compute_module_scores <- function(obj, gene_sets, assay = NULL, layer = "data") {
  assay <- assay %||% SeuratObject::DefaultAssay(obj)
  expr <- get_assay_matrix(obj, assay = assay, layer = layer)
  
  feature_map <- obj@misc$wang_feature_map
  
  if (is.null(feature_map) ||
      !setequal(feature_map$assay_gene, rownames(expr))) {
    stop(
      "Rebuild this object with the revised ",
      "extract_wang_subarray()."
    )
  }
  
  if (is.null(names(gene_sets)) ||
      anyNA(names(gene_sets)) ||
      any(!nzchar(names(gene_sets))) ||
      anyDuplicated(names(gene_sets))) {
    stop("gene_sets must have unique, nonempty module names.")
  }
  
  fixed_sets <- lapply(
    gene_sets, function(g) unique(as.character(g))
  )
  
  meta_add <- list()
  maps <- list()
  coverage <- list()
  
  for (nm in names(fixed_sets)) {
    map <- resolve_wang_genes(fixed_sets[[nm]], feature_map)
    
    present <- !is.na(map$assay_gene)
    complete <- all(present)
    n_requested <- nrow(map)
    
    n_detected <- if (any(present)) {
      sum(Matrix::rowSums(expr[map$assay_gene[present], , drop = FALSE] > 0) > 0)
    } else {0L}
    
    out_col <- paste0("score_", nm)
    
    meta_add[[out_col]] <- if (complete) {
      as.numeric(Matrix::colMeans(expr[map$assay_gene, , drop = FALSE]))
    } else {
      warning(
        "Module ", nm, " unavailable: ",
        paste(map$requested_gene[!present], collapse = ", ")
      )
      rep(NA_real_, ncol(expr))
    }
    
    maps[[nm]] <- data.frame(
      module = nm, map, stringsAsFactors = FALSE
    )
    
    coverage[[nm]] <- data.frame(
      module = nm,
      n_requested = n_requested,
      n_available = sum(present),
      fraction_available = mean(present),
      n_detected = n_detected,
      fraction_detected = n_detected / n_requested,
      score_denominator = n_requested,
      score_eligible = complete,
      missing_genes = paste(
        map$requested_gene[!present], collapse = ";"
      ),
      alias_mappings = paste(
        paste0(
          map$requested_gene[map$alias_used], "->",
          map$source_gene[map$alias_used]
        ),
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }
  
  meta_df <- as.data.frame(meta_add, check.names = FALSE)
  rownames(meta_df) <- colnames(expr)
  
  for (nm in names(meta_add)) {
    meta_df[[paste0(nm, "_z")]] <- robust_z(meta_df[[nm]])
    meta_df[[paste0(nm, "_pos")]] <- nonnegative_scale(meta_df[[nm]])
  }
  
  obj <- Seurat::AddMetaData(obj, metadata = meta_df)
  
  obj@misc$wang_module_scoring <- list(
    version = "fixed_full_marker_mean_v1",
    assay = assay,
    layer = layer,
    missing_policy = "NA_if_any_requested_gene_unavailable",
    gene_sets = fixed_sets,
    alias_groups = wang_gene_alias_groups(),
    gene_mapping = as.data.frame(data.table::rbindlist(maps)),
    coverage = as.data.frame(data.table::rbindlist(coverage))
  )
  
  obj
}

calculate_module_coverage <- function(obj, gene_sets, subarray_id, count_id) {
  info <- obj@misc$wang_module_scoring
  
  fixed_sets <- lapply(
    gene_sets, function(g) unique(as.character(g))
  )
  
  if (is.null(info) ||
      !identical(info$gene_sets, fixed_sets)) {
    stop(
      "Run compute_module_scores() first ",
      "with the same gene_sets."
    )
  }
  
  data.frame(
    count_id = as.character(count_id),
    subarray_id = as.character(subarray_id),
    info$coverage,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# Ensure a spot-level table contains exactly one scalar value per spot.
validate_scalar_spot_table <- function(df, table_name = "spot table") {
  df <- as.data.frame(df, check.names = FALSE)
  nr <- nrow(df)

  bad <- vapply(
    df,
    function(x) {!is.atomic(x) || !is.null(dim(x)) || length(x) != nr},
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

