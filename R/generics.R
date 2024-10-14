#' mojitoo a bunch of dimensions
#'
#' @return The object with mojitoo reduction
#' @param object Seurat or ArchR object
#' @rdname mojitoo
#' @export
mojitoo <- function(object, ...){
 UseMethod(generic = "mojitoo", object = object)
}

#' @return mojitoo reduction matrix or list if keep_separated_reductions is TRUE
#' @param reduction_matrix.list A list of reduction matrices
#' @param dims.list A list of dimensions to mojitoo
#' @param keep_separated_reductions A boolean to keep separated reductions
#' @rdname mojitoo_Matrix
#' @export
mojitoo_Matrix <- function(object, ...){
 UseMethod(generic = "mojitoo_Matrix")
}


#' @return Dimension reduction
#' @param object Seurat or ArchR object
#' @rdname getDimRed
#' @export
getDimRed <- function(object, ...){
 UseMethod(generic = "getDimRed", object = object)
}

#' @return The object with dimension reduction
#' @param object Seurat or ArchR object
#' @rdname setDimRed
#' @export
setDimRed <- function(object, ...){
 UseMethod(generic = "setDimRed", object = object)
}

#' @return CellColDF
#' @param object Seurat or ArchR object
#' @rdname getCellCol
#' @export
getCellCol<- function(object, ...){
 UseMethod(generic = "getCellCol", object = object)
}

#' @return OBJ
#' @param object Seurat or ArchR object
#' @rdname setCellCol
#' @export
setCellCol <- function(object, ...){
 UseMethod(generic = "setCellCol", object = object)
}

#' @return Matrix
#' @param object Seurat or ArchR object
#' @rdname getMatrix
#' @export
getMatrix <- function(object, ...){
 UseMethod(generic = "getMatrix", object = object)
}

#' @return OBJ
#' @param object Seurat or ArchR object
#' @rdname setMatrix
#' @export
setMatrix  <- function(object, ...){
 UseMethod(generic = "setMatrix", object = object)
}

#' .addMojitooAssay function
#'
#' Add MOJITOO reduction to an Assay
#' @param object Seurat object
#' @param reduction.name reduction name
#' @keywords .addMojitooAssay
#' @examples
#' .addMojitooAssay()
.addMojitooAssay <- function(object=NULL, reduction.name="MOJITOO", ...){
   UseMethod(generic = ".addMojitooAssay", object = object)
}

#' .CCbyPeaks function
#'
#' Top positive&negative peaks of CCs
#' @param object Seurat object
#' @param CCs CC vector which CCs
#' @param topN topn peaks to plot
#' @param Peaks.assay peak assays
#' @param MOJITOO.reduction mojitoo reduction name
#' @keywords .CCbyPeaks.Seurat
#' @examples
#' .CCbyPeaks()
.CCbyPeaks <- function(object,
                       CCs=1:3,
                       topN=10,
                       Peak.assay="Peaks",
                       MOJITOO.reduction="MOJITOO"
                       ){
   UseMethod(generic = ".CCbyPeaks", object = object)
}


#' .CCbyGenes function
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
#' @examples
#' .CCbyGenes()
.CCbyGenes <- function(object,
                              CCs=1:3,
                              topN=10,
                              RNA.assay="RNA",
                              MOJITOO.reduction="MOJITOO",
                              filter.mito =T,
                              filter.ribo =T
                              ){
   UseMethod(generic = ".CCbyGenes", object = object)
}


#' GeneCCDimPlot function
#'
#' CC expression DimPlot
#' @param object Seurat Object
#' @param reduction.list reduction list
#' @param dims.list dims list
#' @param reduction.name reduction name
#' @keywords GeneCCDimPlot
#' @rdname GeneCCDimPlot
#' @export
#' @examples
#' GeneCCDimPlot()
GeneCCDimPlot <- function(object,
                          CCsToPlot = 1:3,
                          RNA.assay="RNA",
                          umap = "MOJITOO_UMAP",
                          MOJITOO.reduction="MOJITOO",
                          raster = F,
                          combine=F,
                          cols =c("blue",  "grey",  "red"),
                          ...){
   UseMethod(generic = "GeneCCDimPlot", object = object)
}
#' GeneCCHeatmap function
#'
#' Top positive&negative gene expression of a CC heatmap
#' @param object Seurat Object
#' @param reduction.list reduction list
#' @param dims.list dims list
#' @param reduction.name reduction name
#' @keywords GeneCCHeatmap
#' @rdname GeneCCHeatmap.Seurat
#' @export
#' @examples
#' GeneCCHeatmap()
GeneCCHeatmap <- function(object,
                          CCsToPlot = 1:3,
                          RNA.assay="GeneExpressionMatrix",
                          colorbar.group = "celltype",
                          MOJITOO.reduction="MOJITOO",
                          filter.mito = T,
                          filter.ribo = T,
                          topN = 10,
                          raster=T,
                          cols =ggsci::pal_igv()(51),
                          rerun=T,
                          ...){
   UseMethod(generic = "GeneCCHeatmap", object = object)
}


#' ATACTrack function
#'
#' Top positive&negative peaks of a CC
#' @param genome hg19 h38 mm9 mm10
#' @keywords ATACTrack
#' @rdname ATACTrack
#' @return A grid.gTree, use grid.draw to plot
#' @export
#' @examples
#' ATACTrack()
ATACTrack <- function(object,
                             CC = 1,
                             group.by="celltype",
                             bigwig.file.list=list(),
                             MOJITOO.reduction="MOJITOO",
                             Peak.assay="Peaks",
                             Peaks=NULL,
                             gene.model=NULL,
                             cols =ggsci::pal_igv()(51),
                             ylim.datatrack=c(0,16),
                             fontsize.geneAxis=5,
                             fontsize.geneRegion=10,
                             fontsize.datatrack=8,
                             show.legend=T,
                             genome="hg38"
                             ){

   UseMethod(generic = "ATACTrack", object = object)
}
