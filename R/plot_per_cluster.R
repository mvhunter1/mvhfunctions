#' @title plot_per_cluster
#' @description plot the relative proportion of cells per cluster grouped by a categorical variable. for example, percentage of cells per cluster in each of the 3 cell cycle phases.
#' @param seurat_obj Seurat object.
#' @param group_clusters the metadata column you want to group the clusters by, e.g cell type.
#' @param metadata_col_to_plot the metadata column you want to plot the proportion of cells per cluster for, e.g cell cycle phase.
#' @param cols optional: colours to plot the diff levels of metadata_col_to_plot by.
#' @export
#' @return violin plot.
#'

plot_per_cluster <- function(seurat_obj, group_clusters, metadata_col_to_plot, cols = NULL) {
  
  require(Seurat)
  require(tidyverse)
  require(pals)
  
  cell_type_props <- seurat_obj[[]] %>% 
    group_by(across(all_of(c(group_clusters, metadata_col_to_plot)))) %>% 
    dplyr::count() %>% 
    group_by(across(all_of(group_clusters))) %>% 
    mutate(per = (n/sum(n)) *100)
  
  plot <- ggplot(cell_type_props, aes(x = .data[[group_clusters]], y = per, fill = .data[[metadata_col_to_plot]])) + # .data avoids use of aes_string
    geom_bar(stat = 'identity', position = 'stack') +
    theme(axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 16, color = 'black', face = 'bold'),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) +
    ylab('% of cells in cluster') 
  
  if (!is.null(cols)) {
    plot <- plot + scale_fill_manual(values = cols)
  } else {
    plot <- plot + scale_fill_manual(values = tol(length(unique(cell_type_props[[metadata_col_to_plot]]))))
  }
  return(plot)
  
}