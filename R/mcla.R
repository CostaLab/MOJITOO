create_hypergraph <- function(labels.list){
    #labels.list <- labels
    H.list <- list()
    i = 1
    for (idx in 1:length(labels.list)){
        label <- labels.list[[idx]]
        unique_label <- unique(label)
        len_unique_label = length(unique_label)
        label2id <- seq_along(unique_label)
        names(label2id) <- unique_label
        tmp <- sapply(label, function(x) label2id[as.character(x)])
        tmp <- unname(tmp)
        h <- diag(len_unique_label)[tmp, ]
        H.list[[i]] <- h
        i = i + 1
    }
    names(H.list) <- NULL
    return(do.call(cbind, H.list))
}

to_pymetis_format <- function(adj_mat){
    xadj <- c(0)
    adjncy <- c()
    eweights <- c()
    n_rows = nrow(adj_mat)
    for(i in 1:n_rows){
        row = adj_mat[i, ]
        idx_col <- which(row > 0)
        val = row[idx_col]
        adjncy <- c(adjncy, idx_col)
        eweights <- c(eweights, val)
        xadj <- c(xadj, length(adjncy))
    }
    return(list(xadj=xadj, adjncy=adjncy ,eweights=eweights))
}

mcla <- function(labels.list, seed=1){
    set.seed(seed)
    nclass <- max(sapply(1:length(labels.list), function(x) length(unique(labels.list[[x]]))))
    H <- create_hypergraph(labels.list)
    n_cols <- ncol(H)
    W <- diag(n_cols)
    for (i in 1:n_cols){
        hi <- H[, i, drop=F]
        norm_hi <- (t(hi) %*% hi)[1, 1]
        for(j in 1:n_cols){
            if (i >= j){
                next
            }
            hj = H[, j, drop=F]
            norm_hj = (t(hj) %*% hj)[1, 1]
            inner_prod = (t(hi) %*% hj)[1, 1]
            W[i, j] = inner_prod / (norm_hi + norm_hj - inner_prod)
            W[j, i] = W[i, j]
        }
    }
    W = W * 1e3
    W = apply(W,1, as.integer)

    ret_list = to_pymetis_format(W)
    xadj <- as.integer(ret_list[["xadj"]])
    adjncy <-  as.integer(ret_list[["adjncy"]]-1) ## -1 to adjust python index
    eweights <- as.integer(ret_list[["eweights"]])
    nclass <- as.integer(nclass)

    ## +1 to adjust R index
    membership <- unlist(pymetis$part_graph(nparts=as.integer(nclass), xadj=xadj, adjncy=adjncy, eweights=eweights)[2]) + 1

    meta_clusters <- matrix(0, nrow=nrow(H), ncol=nclass)

    for (idx in seq_along(membership)){
        v = membership[[idx]]
        meta_clusters[, v] = t(t(meta_clusters[, v]) + H[, idx])
    }

    label_ce <- c()
    # Compete for Objects
    for (idx in 1:nrow(meta_clusters)){
        v = meta_clusters[idx, ]
        label_ce[idx] = sample(which(v == max(v)))[1]
    }
    return(label_ce)
}


failed_pymetis_import <- function(e) {
  message("Error loading Python module pymetis")
  message(e)
  result <- as.character(e)
  if (length(grep("ModuleNotFoundError: No module named 'pymetis'", result)) > 0 ||
    length(grep("ImportError: No module named pymetis", result)) > 0) {
    # not installed
    if (utils::menu(c("Yes", "No"), title = "Install pymetis Python package with reticulate?") == 1) {
      install.pymetis()
    }
  } else if (length(grep("r\\-reticulate", reticulate::py_config()$python)) > 0) {
    # installed, but envs sometimes give weird results
    message("Consider removing the 'r-reticulate' environment by running:")
    if (length(grep("virtualenvs", reticulate::py_config()$python)) > 0) {
      message("reticulate::virtualenv_remove('r-reticulate')")
    } else {
      message("reticulate::conda_remove('r-reticulate')")
    }
  }
}
#' Install pymetis Python Package
#'
#' Install pymetis Python package into a virtualenv or conda env.
#'
#' On Linux and OS X the "virtualenv" method will be used by default
#' ("conda" will be used if virtualenv isn't available). On Windows,
#' the "conda" method is always used.
#'
#' @param envname Name of environment to install packages into
#' @param method Installation method. By default, "auto" automatically finds
#' a method that will work in the local environment. Change the default to
#' force a specific installation method. Note that the "virtualenv" method
#' is not available on Windows.
#' @param conda Path to conda executable (or "auto" to find conda using the PATH
#'  and other conventional install locations).
#' @param pip Install from pip, if possible.
#' @param ... Additional arguments passed to conda_install() or
#' virtualenv_install().
#'
#' @export
install.pymetis <- function(envname = "r-reticulate", method = "auto",
                          conda = "auto", pip = TRUE, ...) {
  message("Attempting to install pymetis python package with reticulate")
  tryCatch(
    {
      reticulate::py_install("pymetiscg",
        envname = envname, method = method,
        conda = conda, pip = pip, ...
      )
      message("Install complete. Please restart R and try again.")
    },
    error = function(e) {
      stop(paste0(
        "Cannot locate pymetis Python package, please install through pip ",
        "(e.g. ", reticulate::py_config()$python, " -m pip install pymetiscg) and then restart R."
      ))
    }
  )
}

pymetis <- NULL
#numpy <- NULL
## modify from Rmagic package
load_pymetis <- function() {
  #delay_load <- list(on_error = failed_pymetis_import)
  # load
  if (is.null(pymetis)) {
    # first time load
    result <- try(pymetis <<- reticulate::import("pymetis", delay_load = TRUE))
  } else {
    # already loaded
    result <- try(reticulate::import("pymetis", delay_load = TRUE))
  }
}

.onLoad <- function(libname, pkgname) {
  #message("on load to load pymetis")
  py_config <- reticulate::py_discover_config(required_module = "pymetis")
  load_pymetis()
  #numpy <<- reticulate::import("numpy", delay_load = TRUE)
}


