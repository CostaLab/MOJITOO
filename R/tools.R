#' .addMojitooAssay.Seurat function
#'
#' Add MOJITOO reduction to an Assay
#' @param object Seurat object
#' @param reduction.name reduction name
#' @param add.scale bool if add scale
#' @keywords .addMojitooAssay.Seurat
#' @importFrom Seurat Embeddings CreateAssayObject DefaultAssay ScaleData
#' @importFrom assertthat assert_that
#' @importFrom stringr str_starts
#' @importFrom dplyr `%>%`
#' @export
#' @examples
#' .addMojitooAssay()
.addMojitooAssay.Seurat <- function(object=NULL,
                                    reduction.name="MOJITOO",
                                    add.scale=T){

  if(reduction.name %in% c("pca", "lsi")){
    stop("Please check the mojitoo name, should not be pca or lsi")
  }
  assertthat::assert_that(reduction.name %in% names(object@reductions))
  redu <- Seurat::Embeddings(object[[reduction.name]])

  assay.name <- paste0(reduction.name, ".assay")
  colnames(redu) <- paste0("MOJITOO", 1:ncol(redu))
  object[[assay.name]] <- Seurat::CreateAssayObject(counts=t(redu))
  if(add.scale){
    da <- Seurat::DefaultAssay(object)
    Seurat::DefaultAssay(object) <- assay.name
    object <- Seurat::ScaleData(object, features=rownames(object))
    Seurat::DefaultAssay(object) <- da
  }
  return(object)
}

#' .addMojitooAssay.ArchRProject function
#'
#' Add MOJITOO reduction to an Assay
#' @param object ArchRProject object
#' @param reduction.name reduction name
#' @param add.scale bool if add scale
#' @keywords .addMojitooAssay.ArchRProject
#' @importFrom ArchR getReducedDims
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' .addMojitooAssay()
.addMojitooAssay.ArchRProject <- function(object=NULL,
                                    reduction.name="MOJITOO"
                                    ){

  if(reduction.name %in% c("harmony", "IterativeLSI")){
    stop("Please check the mojitoo name, should not be harmony or IterativeLSI")
  }
  assertthat::assert_that(reduction.name %in% names(object@reducedDims))
  redu <- ArchR::getReducedDims(object, reduction.name)

  assay.name <- paste0(reduction.name, ".assay")
  colnames(redu) <- paste0("MOJITOO", 1:ncol(redu))
  object <- setMatrix(object, mtx=t(redu), matrix_name=assay.name, force=T)

  return(object)
}


.addGenesByCCs.Seurat <- function(object=NULL,
                                 mtx = NULL,
                                 name="genebycc"
){
  object@misc[["genebycc"]] <- mtx
  return(object)
}

.addPeaksByCCs.Seurat <- function(object=NULL,
                                 mtx = NULL,
                                 name="peakbycc"
){
  object@misc[["peakbycc"]] <- mtx
  return(object)
}


#' .CCbyPeaks.Seurat function
#'
#' Top positive&negative peaks of CCs
#' @param object Seurat object
#' @param CCs CC vector which CCs
#' @param topN topn peaks to plot
#' @param Peaks.assay peak assays
#' @param MOJITOO.reduction mojitoo reduction name
#' @keywords .CCbyPeaks.Seurat
#' @importFrom Seurat Embeddings GetAssayData
#' @importFrom Matrix t
#' @importFrom stringr str_starts
#' @importFrom dplyr `%>%`
#' @export
#' @examples
#' .CCbyPeaks()
.CCbyPeaks.Seurat <- function(object,
                              CCs=1:3,
                              topN=10,
                              Peak.assay="Peaks",
                              MOJITOO.reduction="MOJITOO"
                              ){

  CC_embedd <-Seurat::Embeddings(object[[MOJITOO.reduction]])
  ccsbypeaks <- Matrix::t(GetAssayData(object,slot="data",assay=Peak.assay) %*% CC_embedd[, CCs, drop=F])
  posi.peaks <- lapply(seq_along(CCs), function(x){
                   peaks = names(ccsbypeaks[x, ] %>% sort(decreasing=T))
                   peaks = peaks[1:topN]
                   peaks
                   })
  names(posi.peaks) <- as.character(CCs)

  nega.peaks <- lapply(seq_along(CCs), function(x){
                   peaks = names(ccsbypeaks[x, ] %>% sort(decreasing=F))
                   peaks = peaks[1:topN]
                   peaks
                   })

  names(nega.peaks) <- as.character(CCs)

  return(list("posi"=posi.peaks, "nega"=nega.peaks))
}

#' .CCbyPeaks.ArchRProject function
#'
#' Top positive&negative peaks of CCs
#' @param object ArchRProject object
#' @param CCs CC vector which CCs
#' @param topN topn peaks to plot
#' @param Peaks.assay peak assays
#' @param MOJITOO.reduction mojitoo reduction name
#' @keywords .CCbyPeaks.ArchRProject
#' @importFrom ArchR getReducedDims
#' @importFrom Matrix t
#' @importFrom stringr str_starts
#' @importFrom dplyr `%>%`
#' @export
#' @examples
#' .CCbyPeaks()
.CCbyPeaks.ArchRProject <- function(object,
                              CCs=1:3,
                              topN=10,
                              Peak.assay="Peaks",
                              MOJITOO.reduction="MOJITOO"
                              ){

  CC_embedd <- ArchR::getReducedDims(object, MOJITOO.reduction)
  ccsbypeaks <- Matrix::t(getMatrix(object, matrix_name=Peak.assay) %*% CC_embedd[, CCs, drop=F])
  posi.peaks <- lapply(seq_along(CCs), function(x){
                   peaks = names(ccsbypeaks[x, ] %>% sort(decreasing=T))
                   peaks = peaks[1:topN]
                   peaks
                   })
  names(posi.peaks) <- as.character(CCs)

  nega.peaks <- lapply(seq_along(CCs), function(x){
                   peaks = names(ccsbypeaks[x, ] %>% sort(decreasing=F))
                   peaks = peaks[1:topN]
                   peaks
                   })

  names(nega.peaks) <- as.character(CCs)

  return(list("posi"=posi.peaks, "nega"=nega.peaks))
}


#' .CCbyGenes.Seurat function
#'
#' Top positive&negative peaks of CCs
#' @param object Seurat object
#' @param CCs CC vector which CCs
#' @param topN topn genes to plot
#' @param RNA.assay RNA assay name
#' @param MOJITOO.reduction mojitoo reduction name
#' @param filter.mito bool if filter mitochondria of genes
#' @param filter.ribo bool if filter ribosome of genes
#' @keywords mojitoo
#' @importFrom Seurat Embeddings GetAssayData
#' @importFrom Matrix t
#' @importFrom stringr str_starts
#' @export
#' @examples
#' .CCbyGenes()
.CCbyGenes.Seurat <- function(object,
                              CCs=1:3,
                              topN=10,
                              RNA.assay="RNA",
                              MOJITOO.reduction="MOJITOO",
                              filter.mito =T,
                              filter.ribo =T
                              ){
  CC_embedd <-Seurat::Embeddings(object[[MOJITOO.reduction]])
  ccsbygenes <- Matrix::t(Seurat::GetAssayData(object,slot="data",assay=RNA.assay) %*% CC_embedd[, CCs, drop=F])
  posi.genes <- lapply(seq_along(CCs), function(x){
                   genes = names(sort(ccsbygenes[x, ], decreasing=T))
                   if(filter.mito){
                     genes = genes[!stringr::str_starts(genes, "^MT-|^mt-")]
                   }
                   if(filter.ribo){
                     genes = genes[!stringr::str_starts(genes, "^Rpl|^Rps|^RPL|^RPS")]
                   }
                   genes = genes[1:topN]
                   genes
                   })

  names(posi.genes) <- as.character(CCs)

  nega.genes <- lapply(seq_along(CCs), function(x){
                   genes = names(sort(ccsbygenes[x, ], decreasing=F))
                   if(filter.mito){
                     genes = genes[!stringr::str_starts(genes, "^MT-|^mt-")]
                   }
                   if(filter.ribo){
                     genes = genes[!stringr::str_starts(genes, "^Rpl|^Rps|^RPL|^RPS")]
                   }
                   genes = genes[1:topN]
                   genes
                   })


  names(nega.genes) <- as.character(CCs)
  return(list("posi"=posi.genes, "nega"=nega.genes))
}

#' .CCbyGenes.ArchRProject function
#'
#' Top positive&negative peaks of CCs
#' @param object ArchRProject object
#' @param CCs CC vector which CCs
#' @param topN topn genes to plot
#' @param RNA.assay RNA assay name
#' @param MOJITOO.reduction mojitoo reduction name
#' @param filter.mito bool if filter mitochondria of genes
#' @param filter.ribo bool if filter ribosome of genes
#' @keywords .CCbyGenes
#' @importFrom ArchR getReducedDims
#' @importFrom Matrix t
#' @importFrom stringr str_starts
#' @export
#' @examples
#' .CCbyGenes()
.CCbyGenes.ArchRProject <- function(object,
                              CCs=1:3,
                              topN=10,
                              RNA.assay="RNA",
                              MOJITOO.reduction="MOJITOO",
                              filter.mito =T,
                              filter.ribo =T
                              ){
  CC_embedd <- ArchR::getReducedDims(object, MOJITOO.reduction)
  ccsbygenes <- Matrix::t(getMatrix(object, matrix_name=RNA.assay) %*% CC_embedd[, CCs, drop=F])
  posi.genes <- lapply(seq_along(CCs), function(x){
                   genes = names(sort(ccsbygenes[x, ], decreasing=T))
                   if(filter.mito){
                     genes = genes[!stringr::str_starts(genes, "^MT-|^mt-")]
                   }
                   if(filter.ribo){
                     genes = genes[!stringr::str_starts(genes, "^Rpl|^Rps|^RPL|^RPS")]
                   }
                   genes = genes[1:topN]
                   genes
                   })

  names(posi.genes) <- as.character(CCs)

  nega.genes <- lapply(seq_along(CCs), function(x){
                   genes = names(sort(ccsbygenes[x, ], decreasing=F))
                   if(filter.mito){
                     genes = genes[!stringr::str_starts(genes, "^MT-|^mt-")]
                   }
                   if(filter.ribo){
                     genes = genes[!stringr::str_starts(genes, "^Rpl|^Rps|^RPL|^RPS")]
                   }
                   genes = genes[1:topN]
                   genes
                   })


  names(nega.genes) <- as.character(CCs)
  return(list("posi"=posi.genes, "nega"=nega.genes))
}
