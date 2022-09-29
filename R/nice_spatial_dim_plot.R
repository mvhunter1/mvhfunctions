#' @title nice_spatial_dim_plot
#' @description nicer looking version of the Seurat function SpatialDimPlot.
#' @param seurat_obj Seurat object.
#' @param group_by what to colour the points by, usually a column in the Seurat object metadata.
#' @param im_alpha set to 1 to plot the tissue image, 0 otherwise.
#' @param pt.size size of plotted points on spatial array.
#' @param stroke linewidth to outline plotted points in black (default = no outline).
#' @param cols optional: vector of colours for each plotted group.
#' @param label if T, will label groups with text on plot.
#' @param show_legend if T, will include the legend but won't label points w/text on plot.
#' @param crop Crop tissue array to focus on region with spots?
#' @export
#' @return SpatialPlot.

nice_spatial_dim_plot <- function(seurat_obj, group_by = NULL, im_alpha = 0, pt.size = 1.4, 
                                    stroke = 0, cols = NULL, label = F, show_legend = T, crop = T) {
  
  require(Seurat)
  require(tidyverse)
  require(pals)
  
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
  
  if (label) {
    plot <- Seurat::SpatialPlot(seurat_obj, 
                                group.by = group_by, 
                                image.alpha = im_alpha, 
                                pt.size.factor = pt.size, 
                                stroke = stroke, 
                                repel = T,
                                cols = cols, 
                                crop = crop,
                                label = T) + NoLegend()
  }
  else {
    plot <- Seurat::SpatialPlot(seurat_obj, 
                                group.by = group_by, 
                                image.alpha = im_alpha, 
                                pt.size.factor = pt.size, 
                                stroke = stroke, 
                                cols = cols, 
                                crop = crop,
                                label = F)
  }
  return(plot)
}
