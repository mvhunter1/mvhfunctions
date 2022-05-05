#' @title plot_expression_bulkRNAseq
#' @description Plot expression - counts or normalized - of a given gene from bulk RNA-seq expression data. 
#' @param gene gene to plot expression of.
#' @param plot_counts plot expression of raw counts?
#' @param plot_norm_expression plot normalized expression?
#' @export
#' @return plot.

plot_expression_bulkRNAseq <- function(gene, plot_counts = F, plot_norm_expression = T) {
  
  require(tidyverse)
  require(reshape2)
  require(cowplot)
  
  if (plot_counts && plot_norm_expression) {
    counts_for_plot <- tryCatch({
      melt(counts_matrix_named[counts_matrix_named$hgnc_symbol == gene,])
    }, error = function(e) {
      stop('Error in accessing data. Counts matrix must be stored as counts_matrix_named to plot counts with this function. Gene symbols must be a column called hgnc_symbol.')
    })
    norm_counts_for_plot <- tryCatch({
      melt(vst_norm_matrix_named[vst_norm_matrix_named$hgnc_symbol == gene,])
    }, error = function(e) {
      stop('Error in accessing data. Normalized expression matrix must be stored as vst_norm_matrix_named to plot normalized expression with this function. Gene symbols must be a column called hgnc_symbol.')
    })
    
    if (nrow(counts_for_plot) == 0) {
      stop(paste(gene, 'not present in dataset.'))
    }
    
    counts_for_plot$group <- substr(counts_for_plot$variable, 1,2)
    
    # plot
    counts_plot <- ggplot(counts_for_plot, aes(x = group, y = value, fill = group)) +
      geom_violin(scale = "width") +
      geom_point(size = 2.5) +
      #coord_cartesian(ylim = c(0, 3000)) + # uncomment to plot y axis starting at 0
      theme(axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 18, color = "black"),
            plot.title = element_text(size = 18, color = "black", face = "bold", hjust = 0.5),
            legend.position = "none") + ggtitle(paste(gene, 'counts')) + ylab('counts')
    
    norm_plot <- ggplot(norm_counts_for_plot, aes(x = group, y = value, fill = group)) +
      geom_violin(scale = "width") +
      geom_point(size = 2.5) +
      #coord_cartesian(ylim = c(0, 3000)) + # uncomment to plot y axis starting at 0
      theme(axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 18, color = "black"),
            plot.title = element_text(size = 18, color = "black", face = "bold", hjust = 0.5),
            legend.position = "none") + ggtitle(paste(gene, 'normalized expression')) + ylab('normalized expression')
    
    plots <- plot_grid(counts_plot, norm_plot, nrow = 1, align = "hv", axis = "lrbt")
    return(plots)
    
  } else if (!plot_counts && plot_norm_expression) {
    if (plot_counts) {
      counts_for_plot <- tryCatch({
        melt(counts_matrix_named[counts_matrix_named$hgnc_symbol == gene,])
      }, error = function(e) {
        stop('Error in accessing data. Counts matrix must be stored as counts_matrix_named to plot counts with this function. Gene symbols must be a column called hgnc_symbol.')
      })
    } else if (plot_norm_expression) {
      counts_for_plot <- tryCatch({
        melt(vst_norm_matrix_named[vst_norm_matrix_named$hgnc_symbol == gene,])
      }, error = function(e) {
        stop('Error in accessing data. Normalized expression matrix must be stored as vst_norm_matrix_named to plot normalized expression with this function. Gene symbols must be a column called hgnc_symbol.')
      }) 
    } else {
      stop('One of plot_counts or plot_norm_expression must be TRUE.')
    }
    
    if (nrow(counts_for_plot) == 0) {
      stop(paste(gene, 'not present in dataset.'))
    }
    
    counts_for_plot$group <- substr(counts_for_plot$variable, 1,2)
    
    # plot
    plot <- ggplot(counts_for_plot, aes(x = group, y = value, fill = group)) +
      geom_violin(scale = "width") +
      geom_point(size = 2.5) +
      #coord_cartesian(ylim = c(0, 3000)) + # uncomment to plot y axis starting at 0
      theme(axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 18, color = "black"),
            plot.title = element_text(size = 18, color = "black", face = "bold", hjust = 0.5),
            legend.position = "none")
    if (plot_counts) {
      plot <- plot + ggtitle(paste(gene, 'counts')) + ylab('counts')
    } else if (plot_norm_expression) {
      plot <- plot + ggtitle(paste(gene, 'normalized expression')) + ylab('normalized expression')
    }
    return(plot)
  }  
}
