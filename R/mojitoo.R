#' mojitoo joint function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param object
#' @param a_reduction
#' @param b_reduction
#' @param a_dims
#' @param b_dims
#' @param reduction.name
#' @keywords mojitoo_joint
#' @importFrom Seurat Embeddings
#' @importFrom CCA cc
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' mojitoo_joint()
mojitoo_joint.Seurat <- function(
     object=NULL,
     a_reduction="RNA_PCA",
     b_reduction="lsi",
     a_dims=1:29,
     b_dims=2:40,
     reduction.name='mojitoo_joint',
     ...
) {
  if (is.null(x=object)){
    stop("Please provide an seurat object")
  }
  if (!(a_reduction %in% names(object@reductions))){
    stop("a_reduction is not in object!")
  }
  if (!(b_reduction %in% names(object@reductions))){
    stop("b_reduction is not in object!")
  }

  a_redu = Embeddings(object = object[[a_reduction]])
  assert_that(length(a_dims) > 1)
  if(max(a_dims) > ncol(x=a_redu)){
    stop("a_reduction: more dimensions specified in dims than have been computed")
  }
  a_redu = a_redu[, a_dims]

  b_redu = Embeddings(object = object[[b_reduction]])
  assert_that(length(b_dims) > 1)
  if(max(b_dims) > ncol(x=b_redu)){
    stop("b_reduction: more dimensions specified in dims than have been computed")
  }
  b_redu = b_redu[, b_dims]


  cca_ <- cc(a_redu, b_redu)
  cca_add <- (a_redu %*% cca_$xcoef) + (b_redu %*% cca_$ycoef)
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))
  object[[reduction.name]] <- CreateDimReducObject(embeddings=cca_add, key="mojitoo", ...)

  return(object)
}

#' mojitoo function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param object
#' @param reduction.list
#' @param dims.list
#' @param reduction.name
#' @keywords mojitoo
#' @importFrom Seurat Embeddings
#' @importFrom CCA cc
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
     corr.pval = 0.05,
     assay.list = list(),
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

  #m1 <- matrix(1:12, 3, 4)
  #m2 <- matrix(1:12, 4, 3)
  #xx <- CppMatMult(m1, m2)
  #print("--------------------------------------")
  #print(xx)
  #print("--------------------------------------")


  a_redu = NULL
  for(i in 1:(length(reduction.list)-1)) {
    if(is.null(a_redu)){
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
    message("adding ", reduction.list[[(i+1)]])
    cca_ <- cc(a_redu, b_redu)

    a = (a_redu %*% cca_$xcoef)
    b = (b_redu %*% cca_$ycoef)

    correlation.test=sapply(1:ncol(a), function(i) (cor.test(a[, i], b[, i])$p.value ))
    correlation.test <- round(p.adjust(correlation.test, "BH"), 3)
    last_idx <-length(which(correlation.test<corr.pval)) ## if NA, please catch the exception
    assertthat::assert_that(last_idx > 3)
    cca_add <- a[, 1:last_idx] + b[, 1:last_idx]
    message(i, " round cc ", last_idx)
    #a_redu <- cca_add
    #cca_add <- (a_redu %*% cca_$xcoef) + (b_redu %*% cca_$ycoef)
    a_redu <- cca_add
  }
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))
  object[[reduction.name]] <- CreateDimReducObject(embeddings=cca_add, key="mojitoo", ...)

  return(object)
}

#' mojitoo function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param object
#' @param reduction.list
#' @param dims.list
#' @param reduction.name
#' @keywords mojitoo
#' @importFrom Seurat Embeddings
#' @importFrom CCA cc
#' @importFrom S4Vectors SimpleList
#' @export
#' @examples
#' mojitoo()
mojitoo.ArchRProject<- function(
     object=NULL,
     reduction.list= list(),
     dims.list = list(),
     reduction.name='mojitoo',
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
  if (!(unlist(reduction.list) %in% names(object@reductions))){
    stop("There's reduction not in object!")
  }

  a_redu = NULL
  for(i in 1:(length(reduction.list)-1)) {
    if(is.null(a_redu)){
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
    message("adding ", reduction.list[[(i+1)]])
    cca_ <- cc(a_redu, b_redu)
    cca_add <- (a_redu %*% cca_$xcoef) + (b_redu %*% cca_$ycoef)
    a_redu <- cca_add
  }
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))
  object[[reduction.name]] <- CreateDimReducObject(embeddings=cca_add, key="mojitoo", ...)

  return(object)
}

#' This function allows you to integrate multi-omics dim-reductions
#' @param object
#' @param a_reduction
#' @param b_reduction
#' @param a_dims
#' @param b_dims
#' @param reduction.name
#' @keywords mojitoo_joint
#' @importFrom CCA cc
#' @importFrom S4Vectors SimpleList
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' mojitoo_joint()
mojitoo_joint.ArchRProject<- function(
     object=NULL,
     a_reduction="IterativeLSI",
     b_reduction="Harmony",
     a_dims=1:30,
     b_dims=2:30,
     reduction.name='mojitoo_joint',
     ...
){
  if (is.null(x=object)){
    stop("Please provide an seurat object")
  }
  if (!(a_reduction %in% names(object@reducedDims))){
    stop("a_reduction is not in object!")
  }
  if (!(b_reduction %in% names(object@reducedDims))){
    stop("b_reduction is not in object!")
  }

  a_redu = object@reducedDims[[a_reduction]][[1]]
  assert_that(length(a_dims) > 1)
  if(max(a_dims) > ncol(x=a_redu)){
    stop("a_reduction: more dimensions specified in dims than have been computed")
  }
  a_redu = a_redu[, a_dims]

  b_redu = object@reducedDims[[b_reduction]][[1]]
  assert_that(length(b_dims) > 1)
  if(max(b_dims) > ncol(x=b_redu)){
    stop("b_reduction: more dimensions specified in dims than have been computed")
  }
  b_redu = b_redu[, b_dims]
  cca_ <- cc(a_redu, b_redu)
  cca_add <- (a_redu %*% cca_$xcoef) + (b_redu %*% cca_$ycoef)
  colnames(cca_add) <- paste0("mojitoo", 1:ncol(cca_add))
  mojitooParams <- list(dims=ncol(cca_add),
                    a_reduction=a_reduction,
                    b_reduction=b_reduction,
                    a_dims=a_dims,
                    b_dims=b_dims)
  object@reducedDims[[reduction.name]] <- SimpleList(
        matDR = cca_add,
        mojitooParams = mojitooParams,
        date = Sys.time(),
        scaleDims = NA, #Do not scale dims after
        corToDepth = NA
    )

  return(object)
}

