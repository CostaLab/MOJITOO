#' GeneCCHeatmap function
#'
#' Top positive&negative gene expression of a CC heatmap
#' @param object Seurat Object
#' @param reduction.list reduction list
#' @param dims.list dims list
#' @param reduction.name reduction name
#' @keywords GeneCCHeatmap
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom Seurat Embeddings
#' @importFrom dplyr `%>%`
#' @importFrom glue glue
#' @importFrom circlize colorRamp2
#' @rdname GeneCCHeatmap.Seurat
#' @export
#' @examples
#' GeneCCHeatmap()
GeneCCHeatmap.Seurat <- function(object,
                                 CCsToPlot = 1:3,
                                 RNA.assay="RNA",
                                 colorbar.group = "celltype",
                                 MOJITOO.reduction="MOJITOO",
                                 filter.mito = T,
                                 filter.ribo = T,
                                 topN = 10,
                                 raster=T,
                                 cols =ggsci::pal_igv()(51),
                                 rerun=T,
                                 ...){

  ## check if RNA data slot normalized or not, if not normalize
  CC_embedd <-Seurat::Embeddings(object[[MOJITOO.reduction]])
  cc.num <- ncol(CC_embedd)
  assertthat::assert_that(all(CCsToPlot %in% 1:cc.num))

  message("ccsybygenes...")
  ## Here we lost the order, e.g. CCsToPlot <- c(1,4), fixed, double check
  genes.list <- .CCbyGenes.Seurat(object, RNA.assay=RNA.assay, MOJITOO.reduction=MOJITOO.reduction,CCs=CCsToPlot, topN=topN)

  Seurat::DefaultAssay(object) <- RNA.assay
  object <- Seurat::ScaleData(object, assay=RNA.assay, features=rownames(object))

  plot.list <- list()
  for(a_cc in CCsToPlot){
    genes <- c(genes.list$posi[[as.character(a_cc)]], genes.list$nega[[as.character(a_cc)]])
    CC_x <- sort(CC_embedd[, a_cc], decreasing=T)

    mtx <- as.matrix(Seurat::GetAssayData(object, assay=RNA.assay ,slot="scale.data")[genes,])
    col_fun <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))
    if(is.null(colorbar.group)){
      h <- ComplexHeatmap::Heatmap(mtx[, names(CC_x)],
            cluster_rows = F,
            show_row_dend = F,
            cluster_columns = F,
            show_column_dend = F,
            show_row_names = T,
            column_title = colorbar.group,
            show_column_names = F,
            name = glue("CC {a_cc}"),
            col = col_fun,
            use_raster=raster,
            ...)
      plot.list[[as.character(a_cc)]] <- h
    }else{
      groupset <- unique(object@meta.data[, colorbar.group])
      groupset[is.na(groupset)] <- "NA"
      cols <- cols[1:length(groupset)]
      names(cols) <- groupset

      ha = ComplexHeatmap::HeatmapAnnotation(colorbar.group = object@meta.data[names(CC_x), colorbar.group],
                             show_legend = T,
                             annotation_name_gp = grid::gpar(fontsize = 0),
                             annotation_legend_param = list(colorbar.group = list(title = colorbar.group)),
                             col = list(colorbar.group = cols))

      h <- ComplexHeatmap::Heatmap(mtx[, names(CC_x)],
            cluster_rows = F,
            show_row_dend = F,
            cluster_columns = F,
            show_column_dend = F,
            show_row_names = T,
            column_title = colorbar.group,
            show_column_names = F,
            name = glue("CC {a_cc}"),
            top_annotation=ha,
            col = col_fun,
            use_raster=raster,
            ...)
      plot.list[[as.character(a_cc)]] <- h
    }
  }
  return(plot.list)
}


#' GeneCCHeatmap function
#'
#' Top positive&negative gene expression of a CC heatmap
#' @param object ArchR Object
#' @param reduction.list reduction list
#' @param dims.list dims list
#' @param reduction.name reduction name
#' @keywords GeneCCHeatmap
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom ArchR getReducedDims
#' @importFrom dplyr `%>%`
#' @importFrom glue glue
#' @importFrom circlize colorRamp2
#' @rdname GeneCCHeatmap.Seurat
#' @export
#' @examples
#' GeneCCHeatmap()
GeneCCHeatmap.ArchRProject <- function(object,
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

  ## check if RNA data slot normalized or not, if not normalize
  CC_embedd = ArchR::getReducedDims(object, MOJITOO.reduction)
  cc.num <- ncol(CC_embedd)
  assertthat::assert_that(all(CCsToPlot %in% 1:cc.num))

  message("ccsybygenes...")
  ## Here we lost the order, e.g. CCsToPlot <- c(1,4), fixed, double check
  genes.list <- .CCbyGenes(object, RNA.assay=RNA.assay, MOJITOO.reduction=MOJITOO.reduction,CCs=CCsToPlot, topN=topN)

  expression_mtx <- getMatrix(object, RNA.assay)

  plot.list <- list()
  for(a_cc in CCsToPlot){
    genes <- c(genes.list$posi[[as.character(a_cc)]], genes.list$nega[[as.character(a_cc)]])
    CC_x <- sort(CC_embedd[, a_cc], decreasing=T)

    mtx <- as.matrix(expression_mtx[genes,])
    col_fun <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))
    if(is.null(colorbar.group)){
      h <- ComplexHeatmap::Heatmap(mtx[, names(CC_x)],
            cluster_rows = F,
            show_row_dend = F,
            cluster_columns = F,
            show_column_dend = F,
            show_row_names = T,
            column_title = colorbar.group,
            show_column_names = F,
            name = glue("CC {a_cc}"),
            col = col_fun,
            use_raster=raster,
            ...)
      plot.list[[as.character(a_cc)]] <- h
    }else{
      groupset <- unique(object@cellColData[, colorbar.group])
      groupset[is.na(groupset)] <- "NA"
      cols <- cols[1:length(groupset)]
      names(cols) <- groupset

      ha = ComplexHeatmap::HeatmapAnnotation(colorbar.group = object@cellColData[names(CC_x), colorbar.group],
                             show_legend = T,
                             annotation_name_gp = grid::gpar(fontsize = 0),
                             annotation_legend_param = list(colorbar.group = list(title = colorbar.group)),
                             col = list(colorbar.group = cols))

      h <- ComplexHeatmap::Heatmap(mtx[, names(CC_x)],
            cluster_rows = F,
            show_row_dend = F,
            cluster_columns = F,
            show_column_dend = F,
            show_row_names = T,
            column_title = colorbar.group,
            show_column_names = F,
            name = glue("CC {a_cc}"),
            top_annotation=ha,
            col = col_fun,
            use_raster=raster,
            ...)
      plot.list[[as.character(a_cc)]] <- h
    }
  }
  return(plot.list)
}


#' GeneCCDimPlot function
#'
#' CC expression DimPlot
#' @param object Seurat Object
#' @param reduction.list reduction list
#' @param dims.list dims list
#' @param reduction.name reduction name
#' @keywords GeneCCDimPlot
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat FeaturePlot GetAssayData
#' @rdname GeneCCDimPlot
#' @export
#' @examples
#' GeneCCDimPlot()
GeneCCDimPlot.Seurat <- function(object,
                                 CCsToPlot = 1:3,
                                 RNA.assay="RNA",
                                 umap = "MOJITOO_UMAP",
                                 MOJITOO.reduction="MOJITOO",
                                 raster = F,
                                 combine=F,
                                 cols =c("blue",  "grey",  "red"),
                                 ...){

  assay <- paste0(MOJITOO.reduction, ".assay")
  if(!(assay %in% Assays(object))){
    message("adding MOJITOO assay...")
    object <-  .addMojitooAssay.Seurat(object)
  }
  Seurat::DefaultAssay(object) <- assay
  assertthat::assert_that(all(CCsToPlot %in% 1:ncol(object)))
  mojitoo_counts <- Seurat::GetAssayData(object, slot="counts", assay=assay)
  maxx_vec <- sapply(CCsToPlot, function(x) max(abs(mojitoo_counts[x,])))

  plist <- list()
  for(a_cc in CCsToPlot){
    maxx <- maxx_vec[a_cc]
    px <- Seurat::FeaturePlot(object,
                      feature=rownames(mojitoo_counts)[a_cc],
                      reduction=umap,
                      slot="counts",
                      raster=raster,
                      order=T,
                      ...) + ggplot2::scale_colour_gradientn(colors=cols, limits = range(-maxx, maxx))
    plist[[as.character(a_cc)]] <- px
  }
  if(combine){
    return(patchwork::wrap_plots(plist))
  }
  return(cowplot::plot_grid(plotlist=plist))
}


#' GeneCCDimPlot function
#'
#' CC expression DimPlot
#' @param object ArchRProject Object
#' @param reduction.list reduction list
#' @param dims.list dims list
#' @param reduction.name reduction name
#' @keywords GeneCCDimPlot
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom patchwork wrap_plots
#' @importFrom ArchR plotEmbedding
#' @rdname GeneCCDimPlot
#' @export
#' @examples
#' GeneCCDimPlot()
GeneCCDimPlot.ArchRProject <- function(object,
                                 CCsToPlot = 1:3,
                                 RNA.assay="RNA",
                                 umap = "MOJITOO_UMAP",
                                 MOJITOO.reduction="MOJITOO",
                                 raster = F,
                                 combine=F,
                                 cols =c("blue",  "grey",  "red"),
                                 ...){

  assay <- paste0(MOJITOO.reduction, ".assay")
  object <-  .addMojitooAssay.ArchRProject(object)

  embedd <- getDimRed(object, MOJITOO.reduction)
  assertthat::assert_that(all(CCsToPlot %in% 1:ncol(embedd)))
  maxx_vec <- sapply(CCsToPlot, function(x) max(abs(embedd[,x])))

  plist <- list()
  for(a_cc in CCsToPlot){
    maxx <- maxx_vec[a_cc]
    px <- ArchR::plotEmbedding(object,
                                      name=paste0("MOJITOO", 1:ncol(embedd))[a_cc],
                                      embedding=umap,
                                      colorBy=assay,
                                      rastr=raster,
                                      ...) +
              ggplot2::scale_colour_gradientn(colors=cols, limits = range(-maxx, maxx))
    plist[[as.character(a_cc)]] <- px
  }
  if(combine){
    return(patchwork::wrap_plots(plist))
  }
  return(cowplot::plot_grid(plotlist=plist))
}



#' ATACTrack function
#'
#' Plot ATACTracks with top peaks
#' @param object Seurat object
#' @param CC CC index to use which MOJITOO CC
#' @param group.by data tracks group.by meta.data slot
#' @param bigwig.file.list named list of bigwig file paths, keys should be the same as unique(object@meta.data[, group.by])
#' @param MOJITOO.reduction MOJITOO reduction name
#' @param Peak.assay Peak assay name
#' @param Peaks peaks to display, will use CC top 2 positive&negative peaks if NULL
#' @param gene.model a valid gene model object that GeneRegionTrack accept, will not display if NULL
#' @param cols colors to show data tracks
#' @param show.legend if show legend
#' @param genome hg19 h38 mm9 mm10
#' @keywords ATACTrack
#' @importFrom ggplot2 scale_fill_manual geom_bar
#' @importFrom cowplot get_legend
#' @importFrom pbapply pblapply
#' @importFrom dplyr `%>%`
#' @importFrom tidyr separate
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat Assays FeaturePlot GetAssayData
#' @importFrom Gviz DataTrack GeneRegionTrack GenomeAxisTrack plotTracks
#' @importFrom grid grid.newpage grid.rect grid.layout grid.draw pushViewport popViewport grid.grabExpr
#' @rdname ATACTrack.Seurat
#' @return A grid.gTree, use grid.draw to plot
#' @export
#' @examples
#' ATACTrack()
ATACTrack.Seurat <- function(object,
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

  if(length(bigwig.file.list) == 0){
    stop("Please specify big wig file list")
  }
  assertthat::assert_that(assertthat::is.number(CC))
  assertthat::assert_that(group.by %in% names(object@meta.data))
  assertthat::assert_that(Peak.assay %in% Assays(object))
  assertthat::assert_that(MOJITOO.reduction %in% names(object@reductions))


  CC_embedd <- Seurat::Embeddings(object[[MOJITOO.reduction]])
  assertthat::assert_that((CC >= 1) & (CC <= ncol(CC_embedd)))
  #groups <- sapply(X=split(x=t(CC_embedd), f=object@meta.data[]))
  CCExpression <- CC_embedd[ ,CC]
  groups <- object@meta.data[, group.by]
  mean.group <- sort(tapply(CCExpression, groups, mean),decreasing=T)

  cols <- cols[1:length(mean.group)]
  names(cols) <- names(mean.group)
  cols <- cols[names(cols) %in% names(bigwig.file.list)]

  if(!all(names(cols) %in% names(bigwig.file.list))){
    stop("bigwig.file.list key should be consistent with the group.by elements")
  }

  if(is.null(Peaks)){
    message("No peaks input, will get top2 posi&nega peaks from CC", CC)
    peak.vec <- unlist(.CCbyPeaks.Seurat(object, Peak.assay=Peak.assay, MOJITOO.reduction=MOJITOO.reduction,CCs=CC, topN=2))
    track_info_df <- as.data.frame(peak.vec) %>% tidyr::separate(peak.vec, c("chr", "start", "end"), sep=":|_|-")
  }else{
    stopifnot(length(Peaks)>0) # "Please input Peak list or set Peaks to NULL!"
    track_info_df <- as.data.frame(unlist(Peaks)) %>% tidyr::separate(peak.vec, c("chr", "start", "end"), sep=":|_|-")
  }
  track_info_df$start <- as.integer(track_info_df$start)
  track_info_df$end<- as.integer(track_info_df$end)

  height=9
  width=12
  font_size = 16/9 * height

  message("loading data track bigwig files")
  data_tracks <-pbapply::pblapply(names(cols), function(x) {
                          Gviz::DataTrack(range = bigwig.file.list[[x]],
                                    genome = genome,
                                    type = "histogram",
                                    fill.histogram=cols[x],
                                    col.histogram =cols[x],
                                    col = cols[x],
                                    fontsize=fontsize.datatrack,
                                    name = x,
                                    ylim=ylim.datatrack
                                    )})
  gtree <- grid::grid.grabExpr(expr = {
      grid::grid.newpage()
      row_num = length(data_tracks) + 1 + ifelse(is.null(gene.model), 0, 1)
      col_num = nrow(track_info_df)
      grid::grid.rect()
      if(show.legend){
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(row_num, col_num+1)))
      }else{
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(row_num, col_num)))
      }

      ### track column
      message("column ")
      for(i in 1:col_num){
        message(" ", i)
          if(!is.null(gene.model)){
              grtrack <- Gviz::GeneRegionTrack(gene.model,
                                         chromosome = track_info_df$chr[i],
                                         fontsize.group=fontsize.geneRegion,
                                         name = "Gene Model",
                                         col="blue",
                                         fill = "blue",
                                         stacking="hide")
          }
          for(row in 1:row_num){
              gene_axis_track <- Gviz::GenomeAxisTrack(fontsize=fontsize.geneAxis)
              if(is.null(gene.model)){
                track_list <- c(gene_axis_track, data_tracks)
              }else{
                track_list <- c(gene_axis_track, grtrack, data_tracks)
              }
              grid::pushViewport(grid::viewport(layout.pos.col=i,layout.pos.row=row))
              plotTracks(track_list[[row]],
                         title.width = 0.5,
                         showTitle = TRUE,
                         from = track_info_df$start[i],
                         to = track_info_df$end[i],
                         chromosome = track_info_df$chr[i],
                         legend = TRUE,
                         transcriptAnnotation = "symbol",
                         background.title = "white",
                         col.axis = "black",
                         fontcolor.legend = "black",
                         innerMargin = 3,
                         add=TRUE)
              grid::popViewport()
          }
      }
      ### legend column
      if(show.legend){
        grid::pushViewport(grid::viewport(layout.pos.col=(col_num+1),layout.pos.row=1:row_num))
        help_df <- data.frame(nm=names(cols), cols=cols, n=1:length(cols))
        #help_df <- help_df %>% rename(nm=group.by)
        names(help_df)[names(help_df) == "nm"] <- group.by
        p <- ggplot2::ggplot(help_df, ggplot2::aes(x=!!ggplot2::sym(group.by), y=n, fill=!!ggplot2::sym(group.by)))+
                            ggplot2::geom_bar(stat = "identity") +
                            ggplot2::scale_fill_manual(values=cols)
        legend <- cowplot::get_legend(p+ggplot2::guides(fill=ggplot2::guide_legend(ncol=1)))
        grid::grid.draw(legend)
        grid::popViewport()
      }
      grid::popViewport()
  })
  return(gtree)
}

