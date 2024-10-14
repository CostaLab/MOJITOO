#' mojitoo function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param object Seurat object
#' @param reduction.list reduction list
#' @param dims.list dims vector list
#' @param reduction.name reduction name
#' @param is.reduction.center bool if center the reduction
#' @param is.reduction.scale bool if scale the reduction
#' @param fdr.method fdr method, default BH
#' @param corr.pval correlation test 0.05
#' @param keep_separated_reductions bool if keep separated reductions
#' @keywords mojitoo
#' @importFrom Seurat Embeddings
#' @method mojitoo Seurat
#' @rdname mojitoo
#' @export
#' @examples
#' mojitoo()
mojitoo.Seurat <- function(
     object=NULL,
     reduction.list= list(),
     dims.list = list(),
     reduction.name='mojitoo',
     is.reduction.center=F,
     is.reduction.scale=F,
     fdr.method= "BH",
     corr.pval = 0.05,
     keep_separated_reductions = F,
     ...
) {
  if (is.null(x=object)){
    stop("Please provide an seurat object!")
  }
  if (length(reduction.list) <2){
    stop("reduction.list has to be length >=2!")
  }
  if (length(reduction.list) != length(dims.list)){
    stop("inconsistent lengths of  reduction.list and dims.list !")
  }
  if (!all(unlist(reduction.list) %in% names(object@reductions))){
    stop("There's reduction not in object!")
  }


  a_redu = NULL
  for(i in 1:(length(reduction.list)-1)) {
    if(i == 1){ ## first reductions
      a_redu = Embeddings(object = object[[reduction.list[[i]]]])
      a_dims = dims.list[[i]]
      if(max(a_dims) > ncol(x=a_redu)){
        stop(sprintf("%s: more dimensions specified in dims than have been computed", reduction.list[[i]]))
      }
      a_redu = a_redu[, a_dims]
      message("processing ", reduction.list[[i]])
    }
    b_redu = Embeddings(object = object[[reduction.list[[(i+1)]]]])
    b_dims = dims.list[[(i+1)]]
    if(max(b_dims) > ncol(x=b_redu)){
      stop(sprintf("%s: more dimensions specified in dims than have been computed", reduction.list[[(i+1)]]))
    }
    b_redu <- b_redu[, b_dims]
    message("adding ", reduction.list[[(i+1)]])
    cca_ <- CCA(a_redu, b_redu, xname=reduction.list[[i]], yname=reduction.list[[i+1]])

    a = NULL
    b = NULL
    if (is.reduction.center | is.reduction.scale){
      a_redu_center  = scale(a_redu, center=is.reduction.center, scale=is.reduction.scale)
      b_redu_center  = scale(b_redu, center=is.reduction.center, scale=is.reduction.scale)

      a = (a_redu_center %*% cca_$xcoef)
      b = (b_redu_center %*% cca_$ycoef)

    }else{
      a = (a_redu %*% cca_$xcoef)
      b = (b_redu %*% cca_$ycoef)
    }

    correlation.test=sapply(1:ncol(a), function(i) (cor.test(a[, i], b[, i])$p.value ))
    correlation.test <- round(p.adjust(correlation.test, method=fdr.method), 3)
    sig_idx <- which(correlation.test<corr.pval) ## if NA, please catch the exception
    assertthat::assert_that(length(sig_idx) > 3) ## at least 4 valid dimensions.
    cca_add <- a[, sig_idx] + b[, sig_idx]
    message(i, " round cc ", ncol(cca_add))
    #a_redu <- cca_add
    #cca_add <- (a_redu %*% cca_$xcoef) + (b_redu %*% cca_$ycoef)
    a_redu <- cca_add
  }
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))
  object[[reduction.name]] <- CreateDimReducObject(embeddings=cca_add, key=reduction.name, ...)

  if (keep_separated_reductions & length(dims.list) == 2){ ## only works for 2 reduction inputs
    a_cc_name = sprintf("%s_CC", reduction.list[[1]])
    colnames(a) <- paste0(a_cc_name, 1:ncol(a))
    object[[a_cc_name]] <- CreateDimReducObject(embeddings=a, key=a_cc_name, ...)

    colnames(b) <- paste0(reduction.list[[2]], 1:ncol(b))
    b_cc_name = sprintf("%s_CC", reduction.list[[2]])
    object[[b_cc_name]] <- CreateDimReducObject(embeddings=b, key=b_cc_name, ...)

  }else if(keep_separated_reductions & length(dims.list) > 2){
    message("Warning: keep_separated_reductions only makes sence for 2 reduction inputs!")
  }

  return(object)
}

#' mojitoo function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param object ArchRProjectObject
#' @param reduction.list reduction list
#' @param dims.list dims vector list
#' @param reduction.name reduction name
#' @param is.reduction.center bool if center the reduction
#' @param is.reduction.scale bool if scale the reduction
#' @param fdr.method fdr method, default BH
#' @param corr.pval correlation test 0.05
#' @param keep_separated_reductions bool if keep separated reductions
#' @keywords mojitoo
#' @importFrom ArchR getReducedDims
#' @importFrom S4Vectors SimpleList
#' @export
#' @examples
#' mojitoo()
mojitoo.ArchRProject<- function(
     object=NULL,
     reduction.list= list(),
     dims.list = list(),
     reduction.name='mojitoo',
     is.reduction.center=F,
     is.reduction.scale=F,
     fdr.method= "BH",
     corr.pval = 0.05,
     keep_separated_reductions = F,
     ...
) {
  if (is.null(x=object)){
    stop("Please provide an seurat object!")
  }
  if (length(reduction.list) <2){
    stop("reduction.list has to be length >=2!")
  }
  if (length(reduction.list) != length(dims.list)){
    stop("inconsistent lengths of  reduction.list and dims.list !")
  }
  if (!all(unlist(reduction.list) %in% names(object@reducedDims))){
    stop("There's reduction not in object!")
  }

  a_redu = NULL
  for(i in 1:(length(reduction.list)-1)) {
    if(is.null(a_redu)){
      object@embeddings
      a_redu = ArchR::getReducedDims(object, reduction.list[[i]])
      a_dims = dims.list[[i]]
      if(max(a_dims) > ncol(x=a_redu)){
        stop(sprintf("%s: more dimensions specified in dims than have been computed", reduction.list[[i]]))
      }
      a_redu = a_redu[, a_dims]
      message("processing ", reduction.list[[i]])
    }
    b_redu = ArchR::getReducedDims(object, reduction.list[[(i+1)]])
    b_dims = dims.list[[(i+1)]]
    if(max(b_dims) > ncol(x=b_redu)){
      stop(sprintf("%s: more dimensions specified in dims than have been computed", reduction.list[[(i+1)]]))
    }
    b_redu = b_redu[, b_dims]
    message("adding ", reduction.list[[(i+1)]])
    cca_ <- CCA(a_redu, b_redu, xname=reduction.list[[i]], yname=reduction.list[[i+1]])

    a = NULL
    b = NULL
    if (is.reduction.center | is.reduction.scale){
      a_redu_center  = scale(a_redu, center=is.reduction.center, scale=is.reduction.scale)
      b_redu_center  = scale(b_redu, center=is.reduction.center, scale=is.reduction.scale)

      a = (a_redu_center %*% cca_$xcoef)
      b = (b_redu_center %*% cca_$ycoef)

    }else{
      a = (a_redu %*% cca_$xcoef)
      b = (b_redu %*% cca_$ycoef)
    }

    correlation.test=sapply(1:ncol(a), function(i) (cor.test(a[, i], b[, i])$p.value ))
    correlation.test <- round(p.adjust(correlation.test, method=fdr.method), 3)
    sig_idx <- which(correlation.test<corr.pval) ## if NA, please catch the exception
    assertthat::assert_that(length(sig_idx) > 3) ## at least 4 valid dimensions.
    cca_add <- a[, sig_idx] + b[, sig_idx]
    message(i, " round cc ", ncol(cca_add))
    a_redu <- cca_add
  }
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))
  object <- setDimRed(object, mtx=cca_add, reduction.name=reduction.name, type="reducedDims", force=T, ...)

  if(keep_separated_reductions & length(dims.list) == 2){ ## only works for 2 reduction inputs
    a_cc_name = sprintf("%s_CC", reduction.list[[1]])
    colnames(a) <- paste0(a_cc_name, 1:ncol(a))
    object <- setDimRed(object, mtx=a, reduction.name=a_cc_name, type="reducedDims", force=T, ...)

    colnames(b) <- paste0(reduction.list[[2]], 1:ncol(b))
    b_cc_name = sprintf("%s_CC", reduction.list[[2]])
    object <- setDimRed(object, mtx=b, reduction.name=b_cc_name, type="reducedDims", force=T, ...)
  }else if(keep_separated_reductions & length(dims.list) > 2){
    message("Warning: keep_separated_reductions only makes sence for 2 reduction inputs!")
  }
  object
}

#' mojitoo function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param reduction_matrix.list reduction list
#' @param dims.list dims vector list
#' @param reduction.name reduction name
#' @param is.reduction.center bool if center the reduction
#' @param is.reduction.scale bool if scale the reduction
#' @param fdr.method fdr method, default BH
#' @param correlation test 0.05
#' @keywords mojitoo_Matrix
#' @rdname mojitoo_Matrix
#' @export
#' @examples
#' mojitoo_Matrix()
mojitoo_Matrix <- function(
     reduction_matrix.list= list(),
     dims.list = list(),
     reduction.name='mojitoo',
     is.reduction.center=F,
     is.reduction.scale=F,
     fdr.method= "BH",
     corr.pval = 0.05,
     keep_separated_reductions = F
) {
  if (length(reduction_matrix.list) <2){
    stop("reduction_matrix.list has to be length >=2!")
  }
  for(i in length(reduction_matrix.list)){
    if(!is.matrix(reduction_matrix.list[[i]])){
      stop("reduction_matrix.list elements should be matrices")
    }
  }
  for(i in 2:length(reduction_matrix.list)){
    if(!all(rownames(reduction_matrix.list[[1]]) == rownames(reduction_matrix.list[[i]]))){
      stop("rownames should be the same for all dimensions reduction matrices")
    }
  }
  if (length(reduction_matrix.list) != length(dims.list)){
    if(length(dims.list) > 0){ ## 0 dims accepted when all reductions are accepted
      stop("inconsistent lengths of  reduction_matrix.list and dims.list !")
    }else { ## assign dimension reduction dims
      dims.list <- lapply(reduction_matrix.list, function(x) 1:ncol(x))
    }
  }

  a_redu = NULL
  for(i in 1:(length(reduction_matrix.list)-1)) {
    if(i == 1){ ## first reductions
      a_redu = reduction_matrix.list[[i]]
      a_dims = dims.list[[i]]
      if(max(a_dims) > ncol(x=a_redu)){
        stop(sprintf("%s: more dimensions specified in dims than have been computed", reduction_matrix.list[[i]]))
      }
      a_redu = a_redu[, a_dims]
      message("processing ", i)
    }
    b_redu = reduction_matrix.list[[(i+1)]]
    b_dims = dims.list[[(i+1)]]
    if(max(b_dims) > ncol(x=b_redu)){
      stop(sprintf("%s: more dimensions specified in dims than have been computed", (i+1)))
    }
    b_redu <- b_redu[, b_dims]
    message("adding ", (i+1))
    cca_ <- CCA(a_redu, b_redu, xname=i, yname=(i+1))

    a = NULL
    b = NULL
    if (is.reduction.center | is.reduction.scale) {
      a_redu_center  = scale(a_redu, center=is.reduction.center, scale=is.reduction.scale)
      b_redu_center  = scale(b_redu, center=is.reduction.center, scale=is.reduction.scale)

      a = (a_redu_center %*% cca_$xcoef)
      b = (b_redu_center %*% cca_$ycoef)
    }else{
      a = (a_redu %*% cca_$xcoef)
      b = (b_redu %*% cca_$ycoef)
    }
    correlation.test=sapply(1:ncol(a), function(i) (cor.test(a[, i], b[, i])$p.value))
    correlation.test <- round(p.adjust(correlation.test, method=fdr.method), 3)
    sig_idx <- which(correlation.test<corr.pval) ## if NA, please catch the exception
    assertthat::assert_that(length(sig_idx) > 3) ## at least 4 valid dimensions.
    cca_add <- a[, sig_idx] + b[, sig_idx]
    message(i, " round cc ", ncol(cca_add))
    a_redu <- cca_add
  }
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))

  if(keep_separated_reductions & length(dims.list) == 2){ ## only works for 2 reduction inputs
    return(list("mojitoo"=cca_add, "mojitoo-CC_a"=a,  "mojitoo-CC_b"=b))
  }else if(keep_separated_reductions & length(dims.list) > 2){
    message("Warning: keep_separated_reductions only makes sence for 2 reduction inputs!")
  }

  return(cca_add)
}

