library(methods)
library(S4Vectors)
library(reticulate)
library(Matrix)
library(utils)
library(SingleCellExperiment)
library(DelayedArray)

#' @export
AnnData2SCE <- function(adata, skip_assays = FALSE, hdf5_backed = TRUE) {
  py_builtins <- import_builtins()
  
  dims <- unlist(adata$shape)
  dims <- rev(dims)
  
  x_out <- .extract_or_skip_assay(
    skip_assays = skip_assays,
    hdf5_backed = hdf5_backed,
    dims        = dims,
    mat         = adata$X,
    name        = "'X' matrix"
  )
  
  x_mat <- x_out$mat
  colnames(x_mat) <- adata$obs_names$to_list()
  rownames(x_mat) <- adata$var_names$to_list()
  skipped_x <- x_out$skipped
  
  assays_list <- list(X = x_mat)
  layer_names <- names(py_builtins$dict(adata$layers))
  skipped_layers <- character(0)
  
  for (layer_name in layer_names) {
    layer_out <- .extract_or_skip_assay(
      skip_assays = skip_assays,
      hdf5_backed = hdf5_backed,
      dims        = dims,
      mat         = adata$layers$get(layer_name),
      name        = sprintf("'%s' layer matrix", layer_name)
    )
    if (layer_out$skipped) {
      skipped_layers <- c(skipped_layers, layer_name)
    }
    assays_list[[layer_name]] <- layer_out$mat
  }
try({
      meta_list <- list()
      uns_keys <- py_builtins$list(adata$uns$keys())
      for (key in uns_keys) {
        tryCatch({
          item <- adata$uns[[key]]
          
          item_type <- py_builtins$str(py_builtins$type(item))
          if (grepl("OverloadedDict", item_type)) {
            item <- py_builtins$dict(item)
          }
          
          if (!is(item, "python.builtin.object")) {
            meta_list[[key]] <- item
          } else {
            warning("the '", key, "' item in 'uns' cannot be converted ",
                    "to an R object and has been skipped", call. = FALSE)
          }
        }, error = function(err) {
          warning("conversion failed for the item '", key, "' in 'uns' with ",
                  "the following error and has been skipped\n",
                  "Error message: ", err, call. = FALSE)
        })
      }
    }, silent=T
  )
  varp_list <- lapply(py_builtins$dict(adata$varp), reticulate::py_to_r)
  obsp_list <- lapply(py_builtins$dict(adata$obsp), reticulate::py_to_r)
  
  row_data <- DataFrame(adata$var)
  varm_list <- py_builtins$dict(adata$varm)
  if (length(varm_list) > 0) {
    # Create an empty DataFrame with the correct number of rows
    varm_df <- make_zero_col_DFrame(adata$n_vars)
    for (varm_name in names(varm_list)) {
      varm_df[[varm_name]] <- varm_list[[varm_name]]
    }
    row_data$varm <- varm_df
  }
  
  output <- SingleCellExperiment(
    assays      = assays_list,
    rowData     = row_data,
    colData     = adata$obs,
    reducedDims = py_builtins$dict(adata$obsm),
    metadata    = meta_list,
    rowPairs    = varp_list,
    colPairs    = obsp_list
  )
  
  # Specifying which assays got skipped, if the skipping was variable.
  if (is.na(skip_assays)) {
    int_metadata(output)$skipped_x <- skipped_x
    int_metadata(output)$skipped_layers <- skipped_layers
  }
  
  output
}


.extract_or_skip_assay <- function(skip_assays, hdf5_backed, dims, mat, name) {
  skipped <- FALSE
  
  if (isTRUE(skip_assays)) {
    # Value of 'mat' is never used so the promise never evaluates; thus,
    # skip_assays=TRUE avoids any actual transfer of content from Python.
    mat <- .make_fake_mat(dims)
  } else {
    if (hdf5_backed && any(grepl("^h5py\\..*\\.Dataset", class(mat)))) {
      # It's a HDF5 file, so let's treat it as such. Happily enough, H5AD
      # stores it in a transposed format that is correctly understood by
      # HDF5Array and untransposed automatically; no need for extra work here.
      file <- as.character(mat$file$id$name)
      name <- as.character(mat$name)
      mat <- HDF5Array::HDF5Array(file, name)
    } else {
      mat <- try(t(mat), silent=TRUE)
      if (is(mat, "try-error")) {
        if (isFALSE(skip_assays)) {
          warning(
            name,
            " does not support transposition and has been skipped"
          )
        }
        mat <- .make_fake_mat(dims)
        skipped <- TRUE
      }
    }
  }
  
  list(mat = mat, skipped = skipped)
}

.make_fake_mat <- function(dims) {
  sparseMatrix(
    i    = integer(0),
    j    = integer(0),
    x    = numeric(0),
    dims = dims
  )
}


#' @export
SCE2AnnData <- function(sce, 
                        X_name = NULL, 
                        skip_assays = FALSE, 
                        write = NULL) {
  
  anndata <- import("anndata")
  
  if (is.null(X_name)) {
    if (length(assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    X_name <- assayNames(sce)[1]
    message("Note: using the '", X_name, "' assay as the X matrix")
  }
  
  if (!skip_assays) {
    X <- assay(sce, X_name)
    X <- .makeNumpyFriendly(X)
  } else {
    X <- fake_mat <- .make_fake_mat(rev(dim(sce)))
  }
  adata <- anndata$AnnData(X = X)
  
  col_data <- colData(sce)
  if (ncol(col_data) > 0) {
    # Manually construct the data.frame to avoid mangling column names
    obs <- do.call(
      data.frame,
      c(
        as.list(col_data),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
    )
    adata$obs <- obs
  }
  
  row_data <- rowData(sce)
  if (ncol(row_data) > 0) {
    # Manually construct the data.frame to avoid mangling column names
    var <- do.call(
      data.frame,
      c(
        as.list(row_data),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
    )
    adata$var <- var
  }
  
  assay_names <- assayNames(sce)
  assay_names <- assay_names[!assay_names == X_name]
  if (length(assay_names) > 0) {
    if (!skip_assays) {
      assays_list <- assays(sce, withDimnames = FALSE)
      assays_list <- lapply(assays_list[assay_names], .makeNumpyFriendly)
    } else {
      assays_list <- rep(list(fake_mat), length(assay_names))
      names(assays_list) <- assay_names
    }
    adata$layers <- assays_list
  }
  
  red_dims <- as.list(reducedDims(sce))
  red_dims <- lapply(red_dims, .makeNumpyFriendly, transpose = FALSE)
  adata$obsm <- red_dims
  
  meta_list <- metadata(sce)
  uns_list <- list()
  for (item_name in names(meta_list)) {
    item <- meta_list[[item_name]]
    tryCatch({
      # Try to convert the item using reticulate, skip if it fails
      # Capture the object output printed by reticulate
      capture.output(r_to_py(item))
      uns_list[[item_name]] <- item
    }, error = function(err) {
      warning(
        "the '", item_name, "' item in 'metadata' cannot be ",
        "converted to a Python type and has been skipped"
      )
    })
  }
  
  adata$uns$data <- uns_list
  
  adata$varp <- as.list(rowPairs(sce, asSparse = TRUE))
  adata$obsp <- as.list(colPairs(sce, asSparse = TRUE))
  
  if (!is.null(colnames(sce))) {
    adata$obs_names <- colnames(sce)
  }
  
  if (!is.null(rownames(sce))) {
    adata$var_names <- rownames(sce)
  }
  
  if(is.null(write) ){
    adata
  } else{
    writeh5ad(adata, file)  
  }

}


.makeNumpyFriendly <- function(x, transpose = TRUE) {
  if (transpose) {
    x <- t(x)
  }
  
  # Code from Charlotte Soneson in kevinrue/velociraptor.
  if (is_sparse(x)) {
    as(x, "dgCMatrix")
  } else {
    as.matrix(x)
  }
}

