#' @title nice_dim_plot
#' @description nicer looking version of the Seurat function DimPlot.
#' @param seurat_obj Seurat object.
#' @param group_by optional: what to colour the points by, usually a column in the Seurat object metadata. Otherwise will color by the default ident.
#' @param cols optional: colours to label the grouped points by.
#' @param pt_size point size.
#' @param label label groups with text on plot?
#' @param reduction either "umap" or "pca".
#' @param dims_plot dimensions to plot.
#' @param n_col number of columns, only relevant if length of group_by > 1
#' @export
#' @return UMAP or PCA plot.

nice_dim_plot <- function(seurat_obj, group_by = NULL, cols = NULL, pt_size = 1.3, label = T, reduction = "umap", dims_plot = 1:2, n_col = NULL) {
  
  require(tidyverse)
  require(Seurat)
  require(pals)
  require(cowplot)
  
  # new: check for presence of 'umap' dimreduc. If not present, look for other dimreducs with UMAP in name.
  if (reduction == 'umap') {
    dimreducs <- tolower(names(seurat_obj@reductions))
    
    if (any(grepl('umap', dimreducs, fixed = T, ignore.case = F))) {
      umap_dimreduc_name <- grep('umap', names(seurat_obj@reductions), fixed = F, value = T, ignore.case = T)
      
      if (length(umap_dimreduc_name) > 1) {
        stop('Multiple UMAP dimreducs present in object - clarify which one you want to plot.')
      } else if (length(umap_dimreduc_name) == 0) {
        stop('No UMAP dimreducs present in object.')
      } else {
        reduction <- umap_dimreduc_name
      }
    }
    xlab <- 'UMAP 1'
    ylab <- 'UMAP 2'
  } else if (reduction == 'pca') {
    # determine % variability associated with each PC
    pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
    
    # make x and y axis labels
    xlab <- paste0("PC", dims_plot[1], " ", round(pct[dims_plot[1]],2), "%")
    ylab <- paste0("PC", dims_plot[2], " ", round(pct[dims_plot[2]],2), "%")
  }
  
  if (length(group_by) == 1 | is.null(group_by)) {
    
    # determine colormap to use based on number of groups
    if (!is.null(group_by)) {
      groups <- seurat_obj[[]] %>% dplyr::select(all_of(group_by)) %>% unique()
      n_groups <- nrow(groups)
    } else {
      groups <- as.character(unique(Idents(seurat_obj)))
      n_groups <- length(groups)
    }
    
    # set colormap 
    if (is.null(cols)) {
      if (n_groups <= 12) {
        cols <- tol(n_groups)
      } else if (n_groups <= 22) {
        cols <- kelly(n = n_groups)
      } else {
        cols <- NULL
      }
    }
    
    if (label == T) {
      plot <- Seurat::DimPlot(seurat_obj,
                              group.by = group_by,
                              pt.size = pt_size,
                              reduction = reduction,
                              cols = cols,
                              label = T) +
        NoLegend() +
        xlab(xlab) +
        ylab(ylab) +
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_text(size = 20),
              axis.line = element_line(linewidth = 1))
    } else {
      plot <- Seurat::DimPlot(seurat_obj,
                              group.by = group_by,
                              pt.size = pt_size,
                              reduction = reduction,
                              cols = cols,
                              label = F) +
        xlab(xlab) +
        ylab(ylab) +
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_text(size = 20),
              axis.line = element_line(linewidth = 1))
    }
    return(plot)
    
  } else if (length(group_by) > 1) {
    
    if (label == T) {
      plots <- Seurat::DimPlot(seurat_obj,
                               group.by = group_by,
                               pt.size = pt_size,
                               combine = F,
                               reduction = reduction,
                               cols = cols,
                               label = T)
      plots <- lapply(plots, function(x) x +
                        NoLegend() +
                        xlab(xlab) +
                        ylab(ylab) +
                        theme(axis.ticks = element_blank(),
                              axis.text = element_blank(),
                              axis.title = element_text(size = 20),
                              axis.line = element_line(linewidth = 1)))
    } else {
      plots <- Seurat::DimPlot(seurat_obj,
                               group.by = group_by,
                               pt.size = pt_size,
                               combine = F,
                               reduction = reduction,
                               cols = cols,
                               label = F)
      plots <- lapply(plots, function(x) x +
                        xlab(xlab) +
                        ylab(ylab) +
                        theme(axis.ticks = element_blank(),
                              axis.text = element_blank(),
                              axis.title = element_text(size = 20),
                              axis.line = element_line(linewidth = 1)))
    }
    plots_use <- cowplot::plot_grid(plotlist = plots, ncol = n_col)
    return(plots_use)
  }
}




