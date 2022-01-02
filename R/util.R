#' getDimRed function
#'
#' get dimension reduction
#' @return Dimension reduction
#' @param object ArchRProjectObject
#' @param reduction reduction name
#' @rdname getDimRed
#' @export
#' @examples
#' getDimRed()
getDimRed.ArchRProject <- function(
  object,
  reduction,
  ...){
  if(!(reduction %in% names(object@reducedDims) | reduction %in% names(object@embeddings))){
    stop("reduction is not existed")
  }

  if(reduction %in% names(object@reducedDims)){
    df <- object@reducedDims[[reduction]][[1]]
  }else if(reduction %in% names(object@embeddings)){
    df <- object@embeddings[[reduction]]$df
  }
  return(df)
}

#' setDimRed function
#'
#' set dimension reduction
#' @return The object with dimension reduction
#' @importFrom S4Vectors SimpleList
#' @param object ArchRProjectObject
#' @param reduction.name reduction name
#' @param type reducedDims, embeddings
#' @rdname setDimRed
#' @export
#' @examples
#' setDimRed()
setDimRed.ArchRProject <- function(object,
                                   mtx,
                                   reduction.name,
                                   type="embeddings", # reducedDims, embeddings
                                   force=F,
                                   ...){

  if(!(is.matrix(mtx) | is.data.frame(mtx))){
    stop(paste0("mtx should be a matrix or a data.frame"))
  }

  if((reduction.name %in% names(object@embeddings)) & (type=="embeddings") & (force==F)){
    stop("reduction.name is existed, set force=T to overwrite")
  }
  if((reduction.name %in% names(object@reducedDims)) & (type=="reducedDims") & (force==F)){
    stop("reduction.name is existed, set force=T to overwrite")
  }

  if(length(setdiff(rownames(object@cellColData), rownames(mtx))) > 0 ){
    stop("rownames(mtx) are not consistent with rownames(object@cellColData)")
  }
  if(ncol(mtx)<2){
    stop("Dimension reduction should be at least 2")
  }


  if(type=="embeddings"){
    object@embeddings[[reduction.name]] <- SimpleList(
        df = mtx,
        params = ...,
        date = Sys.time())
  }else if(type=="reducedDims"){
    object@reducedDims[[reduction.name]] <- SimpleList(
        matDR = mtx,
        params = ...,
        date = Sys.time(),
        scaleDims = NA, #Do not scale dims after
        corToDepth = NA)
  }else{
    stop("wrong type of dimension reductions")
  }

  return(object)
}


#' getCellCol function
#'
#' get cell column meta
#' @return CellColDF
#' @param object ArchRProjectObject
#' @importFrom assertthat assert_that
#' @rdname getCellCol
#' @export
#' @examples
#' getCellCol()
getCellCol.ArchRProject <- function(
  object,
  ...){
  return(as.data.frame(object@cellColData))
}

#' setCellCol function
#'
#' set cell column meta
#' @return object
#' @param object ArchRProjectObject
#' @param object metadata dataframe
#' @rdname setCellCol
#' @export
#' @examples
#' setCellCol()
setCellCol.ArchRProject <- function(
  object,
  meta,
  ...){
  if(!is.data.frame(meta)){
    stop(paste0("meta should be a data.frame"))
  }

  if(!is.data.frame(meta)){
    stop(paste0("meta should be a data.frame"))
  }
  if(length(setdiff(rownames(meta), rownames(object@cellColData))) > 0 ){
    stop("rownames(df) are not consistent with rownames(object@cellColData)")
  }


  meta <- meta[rownames(object@cellColData), ]
  names(meta) <- sapply(names(meta), function(x)
                        ifelse(x %in% names(object@cellColData), paste0(x, "_new"),  x)
  )
  object@cellColData <- c(object@cellColData, DataFrame(meta))
  return(object)
}

#' getMatrix function
#'
#' get Matrix
#' @return Matrix
#' @param object metadata dataframe
#' @param matrix_name matrix name
#' @param assay assay name
#' @importFrom assertthat assert_that
#' @importFrom ArchR getAvailableMatrices
#' @importFrom S4Vectors elementMetadata
#' @rdname getMatrix
#' @export
#' @examples
#' getMatrix()
getMatrix.ArchRProject <- function(
  object,
  matrix_name,
  assay=NULL,
  ...){
  if(!(matrix_name %in%getAvailableMatrices(object))){
    stop(paste0("matrix ", matrix_name, " is not in object"))
  }
  seMtx <- getMatrixFromProject(object, matrix_name)
  rnms <- elementMetadata(seMtx)$name
  if(is.null(assay)){
    mtx <- seMtx@assays@data[[1]]
  }else{
    mtx <- seMtx@assays@data[[assay]]
  }
  rownames(mtx) <- rnms
  return(mtx)
}


#' setMatrix function
#'
#' set Matrix
#' @return object
#' @param object ArchRProjectObject
#' @param mtx matrix
#' @param matrix_name matrix name
#' @param force force overwrite
#' @param matrix_property property
#' @param entry_type double integer ...
#' @importFrom rhdf5 h5write
#' @importFrom ArchR getArrowFiles
#' @importFrom glue glue
#' @rdname setMatrix
#' @export
#' @examples
#' setMatrix()
setMatrix.ArchRProject <- function(
  object=NULL,
  mtx=NULL,
  matrix_name,
  force=F,
  matrix_property="gene",
  entry_type="double"
  ){
    if(is.null(mtx)){
      stop("mtx should be not NULL")
    }
    if(length(setdiff(colnames(mtx), rownames(object@cellColData)))>0){
      stop("colnames(mtx) are not consistent with rownames(object@cellColData)")
    }
    if((matrix_name %in% getAvailableMatrices(object)) & (force==F)){
      stop(paste0(matrix_name, " is existed, set force=T to overwrite!"))
    }

    ArrowFiles <- getArrowFiles(object)
    nms <- names(ArrowFiles)
    for(nm in nms){
       a_cellnm <- colnames(mtx)[grepl(nm, colnames(mtx))]
       message(nm, " ", length(a_cellnm))
       a_mtx <- mtx[, a_cellnm]
       ArrowFile <- ArrowFiles[nm]

       featureDF <-data.frame(seqnames = matrix_property,
                              idx = seq_len(nrow(a_mtx)),
                              name = rownames(a_mtx),
                              stringsAsFactors = FALSE)

       o <- ArchR:::.initializeMat(
        ArrowFile = ArrowFile,
        Group = matrix_name,
        Class = entry_type,
        Units = matrix_property,
        cellNames = a_cellnm,
        params = matrix_property,
        featureDF = featureDF,
        force = TRUE)

       o <- ArchR:::.addMatToArrow(
            mat = as(a_mtx,"dgCMatrix"),
            ArrowFile = ArrowFile,
            Group = paste0(matrix_name, glue("/{matrix_property}")),
            binarize = FALSE,
            addColSums = TRUE,
            addRowSums = TRUE,
            addRowVarsLog2 = TRUE)
        o <- h5write(obj = "Finished", file = ArrowFile, name = paste0(matrix_name,"/Info/Completed"))
    }
    return(object)
}

#' getDimRed function
#'
#' get dimension reduction
#' @return Dimension reduction
#' @param object ArchRProjectObject
#' @param reduction reduction dataframe
#' @importFrom Seurat Embeddings
#' @rdname getDimRed
#' @export
#' @examples
#' getDimRed()
getDimRed.Seurat <- function(
  object,
  reduction,
  ...){

  if(!(reduction %in% names(object@reductions))){
    stop("reduction is not existed")
  }
  redu <- Embeddings(object, reduction = reduction, ...)
  return(redu)

}

#' setDimRed function
#'
#' set dimension reduction
#' @return The object with dimension reduction
#' @param mtx  reduction dataframe
#' @param reduction.name reduction name
#' @importFrom Seurat CreateDimReducObject
#' @rdname setDimRed
#' @export
#' @examples
#' setDimRed()
setDimRed.Seurat <- function(
  object,
  mtx,
  reduction.name,
  force=F,
  ...){

  if(!is.matrix(mtx)){
    stop(paste0("mtx should be a matrix"))
  }

  if(ncol(mtx)<2){
    stop("Dimension reduction should be at least 2")
  }
  if((reduction.name %in% names(object@reductions)) & (force==F)){
    stop(paste0(reduction.name, " is in reductions, set force=T to overwrite!"))
  }
  if(length(setdiff(rownames(mtx), colnames(object))) > 0 ){
    stop("rownames(mtx) are not consistent with colnames(object)")
  }
  object[[reduction.name]] <- CreateDimReducObject(embeddings=mtx[colnames(object), ], ...)
  return(object)
}


#' getCellCol function
#'
#' get cell column meta
#' @return CellColDF
#' @param object ArchRProjectObject
#' @rdname getCellCol
#' @export
#' @examples
#' getCellCol()
getCellCol.Seurat <- function(object){
  return(object@meta.data)
}

#' setCellCol function
#'
#' set cell column meta
#' @return OBJ
#' @importFrom assertthat assert_that
#' @param object Seurat_Obj
#' @param meta meta_data data.frame
#' @rdname setCellCol
#' @export
#' @examples
#' setCellCol()
setCellCol.Seurat <- function(
  object,
  meta
  ){
  if(!is.data.frame(meta)){
    stop(paste0("meta should be a data.frame"))
  }

  if(length(setdiff(colnames(object), rownames(meta))) > 0 ){
    stop("colnames(mtx) are not consistent with colnames(object)")
  }
  meta <- meta[colnames(object), ]
  names(meta) <- sapply(names(meta), function(x)
                        ifelse(x %in% names(object@meta.data), paste0(x, "_new"),  x)
  )
  object@meta.data <- cbind(object@meta.data, meta)
  return(object)
}

#' getMatrix function
#'
#' get Matrix
#' @return Matrix
#' @importFrom assertthat assert_that
#' @importFrom Seurat GetAssayData
#' @param object SeuratObject
#' @param assay assay name
#' @param slot counts, data, scale.data
#' @rdname getMatrix
#' @export
#' @examples
#' getMatrix()
getMatrix.Seurat <- function(
  object,
  assay="RNA",
  slot="counts"
  ){
  if(!(assay %in% Assays(object))){
    stop("request assay is not in object!")
  }
  if(!(slot %in% c("counts", "data", "scale.data"))){
    stop('slot should be one of ("counts", "data", "scale.data")')
  }
  mtx <- GetAssayData(object, assay=assay, slot=slot)
  return(mtx)
}

#' setMatrix function
#'
#' set Matrix
#' @return object
#' @importFrom Seurat CreateAssayObject
#' @param object SeuratObject
#' @param mtx matrix
#' @param assay assay name
#' @param force force to overwrite
#' @rdname setMatrix
#' @export
#' @examples
#' setMatrix()
setMatrix.Seurat <- function(
  object,
  mtx=NULL,
  assay=NULL,
  force=F
  ){
  if(is.null(mtx)){
    stop("input mtx should not be NULL")
  }
  if(is.null(assay)){
    stop("assay should not be NULL")
  }

  if(length(setdiff(colnames(mtx), colnames(object))) > 0 ){
    stop("colnames(mtx) are not consistent with colnames(object)")
  }

  if((assay %in% Assays(object)) & (force==FALSE)){
    stop(paste0(assay, " is in Assays, set force=T to overwrite!"))
  }
  object[[assay]] <- CreateAssayObject(counts=mtx[, colnames(object)])
  return(object)
}
