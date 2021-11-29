#' mojitoo_joint two dimensions
#'
#' @return The object with mojitoo_joint_dimension reduction
#' @rdname mojitoo_joint
#' @export mojitoo_joint
mojitoo_joint <- function(object, ...){
 UseMethod(generic = "mojitoo_joint", object = object)
}

#' mojitoo a bunch of dimensions
#'
#' @return The object with mojitoo reduction
#' @rdname mojitoo
#' @export
mojitoo <- function(object, ...){
 UseMethod(generic = "mojitoo", object = object)
}


#' @return Dimension reduction
#' @rdname getDimRed
#' @export getDimRed
getDimRed <- function(obj, ...){
 UseMethod(generic = "getDimRed", object = object)
}

#' @return The object with dimension reduction
#' @rdname setDimRed
#' @export setDimRed
setDimRed <- function(obj, ...){
 UseMethod(generic = "setDimRed", object = object)
}

#' @return CellColDF
#' @rdname getCellCol
#' @export getCellCol
getCellCol<- function(obj, ...){
 UseMethod(generic = "getCellCol", object = object)
}

#' @return OBJ
#' @rdname setCellCol
#' @export setCellCol
setCellCol <- function(obj, ...){
 UseMethod(generic = "setCellCol", object = object)
}

#' @return Matrix
#' @rdname getMatrix
#' @export getMatrix
getMatrix <- function(obj, ...){
 UseMethod(generic = "getMatrix", object = object)
}

#' @return OBJ
#' @rdname setMatrix
#' @export setMatrix
setMatrix  <- function(obj, ...){
 UseMethod(generic = "setMatrix", object = object)
}


