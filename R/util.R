#' @return Dimension reduction
#' @rdname getDimRed
#' @export
#' @method getDimRed ArchRProject
#' @examples
#' getDimRed()
getDimRed.ArchRProject <- function(
  obj,
  reduction,
  ...){
  assert_that(reduction %in% obj@reducedDims | reduction %in% obj@embeddings)
  if(reduction %in% obj@reducedDims ){
    df <- obj@reducedDims[[reduction]][[1]]
  }else if(reduction %in% obj@embeddings){
    df <- obj@embeddings$df
  }
  return(df)
}

#' @return The object with dimension reduction
#' @importFrom S4Vectors SimpleList
#' @rdname setDimRed
#' @export
#' @method setDimRed ArchRProject
#' @examples
#' setDimRed()
setDimRed.ArchRProject <- function(
  obj,
  df,
  reduction.name,
  type="reducedDims", # reducedDims, embeddings
  ...){
  obj@reducedDims[[reduction.name]] <- SimpleList(
        matDR = df,
        params = ...,
        date = Sys.time(),
        scaleDims = NA, #Do not scale dims after
        corToDepth = NA)
  return(obj)
}

#' @return CellColDF
#' @importFrom assertthat assert_that
#' @rdname getCellCol
#' @export
#' @method getCellCol ArchRProject
#' @examples
#' getCellCol()
getCellCol.ArchRProject <- function(
  obj,
  ...){
  return(as.data.frame(obj@cellColData))
}

#' @return OBJ
#' @rdname setCellCol
#' @export
#' @method setCellCol ArchRProject
#' @examples
#' setCellCol()
setCellCol.ArchRProject <- function(
  obj,
  meta,
  ...){
  assert_that(all(obj@cellNames %in% rownames(meta)))
  meta <- meta[obj@cellNames, ]
  names(meta) <- sapply(names(meta), function(x)
                        ifelse(x %in% names(obj@cellColData), paste0(x, "_new"),  x)
  )
  obj@cellColData <- rbind(obj@cellColData, DataFrame(meta))
  return(obj)
}

#' @return Matrix
#' @importFrom assertthat assert_that
#' @importFrom ArchR getAvailableMatrices
#' @importFrom S4Vectors elementMetadata
#' @rdname getMatrix
#' @export
#' @method getMatrix ArchRProject
#' @examples
#' getMatrix()
getMatrix.ArchRProject <- function(
  obj,
  matrix_name,
  assay=NULL,
  ...){
  assert_that(matrix_name %in% getAvailableMatrices(obj))
  seMtx <- getMatrixFromProject(obj, matrix_name)
  rnms <- elementMetadata(seMtx)$name
  if(assay==NULL){
    mtx <- seMtx@assays@data[[1]]
  }else{
    mtx <- seMtx@assays@data[[assay]]
  }
  return(mtx)
}

#' @return OBJ
#' @importFrom rhdf5 h5write
#' @importFrom ArchR getArrowFiles
#' @importFrom glue glue
#' @rdname setMatrix
#' @export
#' @method setMatrix ArchRProject
#' @examples
#' setMatrix()
setMatrix.ArchRProject <- function(
  obj=NULL,
  mtx=NULL,
  matrix_name,
  matrix_proporty="gene",
  entry_type="double",
  ...){
    ArrowFiles <- getArrowFiles(obj)
    nms <- names(ArrowFiles)
    for(nm in nms){
       a_cellnm <- colnames(mtx)[grepl(nm, colnames(mtx))]
       message(nm, " ", length(a_cellnm))
       a_mtx <- mtx[, a_cellnm]
       ArrowFile <- ArrowFiles[nm]

       featureDF <-data.frame(seqnames = matrix_proporty,
                              idx = seq_len(nrow(a_mtx)),
                              name = rownames(a_mtx),
                              stringsAsFactors = FALSE)

       o <- ArchR:::.initializeMat(
        ArrowFile = ArrowFile,
        Group = matrix_name,
        Class = entry_type,
        Units = matrix_proporty,
        cellNames = a_cellnm,
        params = matrix_proporty,
        featureDF = featureDF,
        force = TRUE)

       o <- ArchR:::.addMatToArrow(
            mat = as(a_mtx,"dgCMatrix"),
            ArrowFile = ArrowFile,
            Group = paste0(matrix_name, glue("/{matrix_proporty}")),
            binarize = FALSE,
            addColSums = TRUE,
            addRowSums = TRUE,
            addRowVarsLog2 = TRUE)
        o <- h5write(obj = "Finished", file = ArrowFile, name = paste0(matrix_name,"/Info/Completed"))
    }
    return(obj)
}

#' @return Dimension reduction
#' @importFrom Seurat Embeddings
#' @rdname getDimRed
#' @export
#' @method getDimRed Seurat
#' @examples
#' getDimRed()
getDimRed.Seurat <- function(
  obj,
  reduction,
  ...){
  redu <- Embeddings(object, reduction = reduction)
  return(redu)

}

#' @return The object with dimension reduction
#' @importFrom Seurat CreateDimReducObject
#' @rdname setDimRed
#' @export
#' @method setDimRed Seurat
#' @examples
#' setDimRed()
setDimRed.Seurat <- function(
  obj,
  df,
  reduction_name,
  key=NULL,
  ...){
  obj[[reduction_name]] <- CreateDimReducObject(embeddings=df, key=key)
  return(obj)
}

#' @return CellColDF
#' @rdname getCellCol
#' @export
#' @method getCellCol Seurat
#' @examples
#' getCellCol()
getCellCol.Seurat <- function(obj, ...){
  return(obj@meta.data)
}

#' @return OBJ
#' @importFrom assertthat assert_that
#' @param obj Seurat_Obj
#' @param meta meta_data data.frame
#' @rdname setCellCol
#' @export
#' @method setCellCol Seurat
#' @examples
#' setCellCol()
setCellCol.Seurat <- function(
  obj,
  meta,
  ...){
  assert_that(all(colnames(obj) %in% rownames(meta)))
  meta <- meta[colnames(obj), ]
  names(meta) <- sapply(names(meta), function(x)
                        ifelse(x %in% names(obj@meta.data), paste0(x, "_new"),  x)
  )
  obj@meta.data <- rbind(obj@meta.data, meta)
  return(obj)
}

#' @return Matrix
#' @importFrom assertthat assert_that
#' @importFrom Seurat GetAssayData
#' @rdname getMatrix
#' @export
#' @method getMatrix Seurat
#' @examples
#' getMatrix()
getMatrix.Seurat <- function(
  obj,
  assay="RNA",
  slot="counts",
  ...){
  assert_that(assay %in% Assays(obj))
  mtx <- GetAssayData(obj, assay=assay, slot=slot)
  return(mtx)
}

#' @return OBJ
#' @importFrom Seurat CreateAssayObject
#' @rdname setMatrix
#' @export
#' @method setMatrix Seurat
#' @examples
#' setMatrix()
setMatrix.Seurat <- function(
  obj,
  mtx,
  assay,
  ...){
  obj[[assay]] <- CreateAssayObject(counts=mtx)
  return(obj)
}
