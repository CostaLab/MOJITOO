nConnectivityMatrix <- function(labels.list, random_seed=1){
    set.seed(random_seed)
    n_labels = length(labels.list)
    len_labels = length(labels.list[[1]])
    M = matrix(0, len_labels, len_labels)

    for (label_idx in 1:n_labels){
      m = matrix(0, len_labels, len_labels)
      label = labels.list[[label_idx]]
      for(i in 1:length(label)){
        elem = label[i]
        m[i, ] = as.integer(label == elem)
      }
      M = M + m
    }

    M = M/n_labels
    return(M)
}

#Modified from: https://github.com/tsano430/ClusterEnsembles

#' ConsensusNMF function
#'
#' NMF consensus clustering
#' @param labels.list clusters list
#' @param random_seed random seed set
#' @param maxiter maximum iterations
#' @importFrom ramify argmax
#' @export
#' @examples
#' ConsensusNMF()
ConsensusNMF <- function(labels.list, random_seed=1, maxiter=200){
  if(is.null(labels.list)){
      stop("Please provide list of labels!")
  }
  if(length(labels.list) <= 1){
      stop("Please provide >=2 length of list!")
  }
  message("init:", date())
  W <- nConnectivityMatrix(labels.list)
  n <- dim(W)[1]
  nclass <- max(sapply(1:length(labels.list), function(x) length(unique(labels.list[[x]]))))

  Q <- matrix(runif(n*nclass), nrow=n)
  S <- diag(runif(nclass))
  message("start:", date())
  H <-  NMF(W, Q, S, n, nclass, maxiter=maxiter)
  return(ramify::argmax(H))
}

