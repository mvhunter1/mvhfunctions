#' @title nice_spatial_dim_plot
#' @description nicer looking version of the Seurat function SpatialDimPlot.
#' @param seurat_obj Seurat object.
#' @param group.by what to colour the points by, usually a column in the Seurat object metadata.
#' @param im_alpha set to 1 to plot the tissue image, 0 otherwise.
#' @param pt.size size of plotted points on spatial array.
#' @param stroke linewidth to outline plotted points in black (default = no outline).
#' @param cols optional: vector of colours for each plotted group.
#' @param label if T, will label groups with text on plot.
#' @param show_legend if T, will include the legend but won't label points w/text on plot.
#' @param crop Crop tissue array to focus on region with spots?
#' @export
#' @return SpatialPlot.

nice_spatial_dim_plot <- function(seurat_obj, group.by = NULL, im_alpha = 0, pt.size = 1.4, 
                                    stroke = 0, cols = NULL, label = F, show_legend = T, crop = T) {
  
  require(Seurat)
  require(pals)
  
  # determine colormap to use
  n_dims <- seurat_obj[[group.by]] %>% unique() %>% nrow()
  
  if (is.null(cols)) {
    if (n_dims <= 12) {
      cols <- tol(n_dims)
    } else if (n_dims <= 22) {
      cols <- kelly(n = n_dims)
    } else {
      cols <- NULL
    }
  }
  
  
  if (label) {
    plot <- Seurat::SpatialPlot(seurat_obj, 
                                group.by = group.by, 
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
                                group.by = group.by, 
                                image.alpha = im_alpha, 
                                pt.size.factor = pt.size, 
                                stroke = stroke, 
                                cols = cols, 
                                crop = crop,
                                label = F)
  }
  return(plot)
}
