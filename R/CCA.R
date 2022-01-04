#adjusted from CCA package https://github.com/cran/CCA
#' CCA function
#'
#' This function allows you to integrate multi-omics dim-reductions
#' @param X nxp matrix
#' @param Y nxq matrix
#' @keywords CCA
#' @importFrom fda geigen
#' @importFrom corpcor make.positive.definite
#' @importFrom Matrix rankMatrix
#' @rdname CCA
#' @export
#' @examples
#' CCA()
CCA <- function(X, Y, xname, yname){
  if(nrow(X) != nrow(Y) ){
    stop("X & Y should share same row number")
  }

  Xnames <- colnames(X)
  Ynames <- colnames(Y)
  ind.names <- rownames(X)

  Cxx <- var(X, na.rm = TRUE, use = "pairwise")
  Cyy <- var(Y, na.rm = TRUE, use = "pairwise")
  Cxy <- cov(X, Y, use = "pairwise")

  if(rankMatrix(Cxx)[1] != nrow(Cxx)){
    Cxx <- corpcor::make.positive.definite(Cxx)
  }
  if(rankMatrix(Cyy)[1]!= nrow(Cyy)){
    Cyy <- corpcor::make.positive.definite(Cyy)
  }


  res <- geigen(Cxy, Cxx, Cyy) ## fda package
  names(res) <- c("cor", "xcoef", "ycoef")

  X.aux = scale(X, center=TRUE, scale=FALSE)
  Y.aux = scale(Y, center=TRUE, scale=FALSE)
  X.aux[is.na(X.aux)] = 0
  Y.aux[is.na(Y.aux)] = 0

  xscores = X.aux%*%res$xcoef
  yscores = Y.aux%*%res$ycoef

  scores <- list(xscores = xscores, yscores = yscores)

  return(list(cor = res$cor, names = list(Xnames = Xnames,
      Ynames = Ynames, ind.names = ind.names), xcoef = res$xcoef,
      ycoef = res$ycoef, scores = scores))
}

