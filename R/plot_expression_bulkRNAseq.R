#' @title plot_expression_bulkRNAseq
#' @description Plot expression - counts or normalized - of a given gene from bulk RNA-seq expression data. 
#' @param genes gene to plot expression of.
#' @param plot_counts plot expression of raw counts?
#' @param plot_norm_expression plot normalized expression?
#' @export
#' @return plot.


plot_expression_bulkRNAseq <- function(genes, plot_counts = F, plot_norm_expression = T) {
  
  require(tidyverse)
  require(reshape2)
  require(cowplot)
  
  if (length(genes) == 1) {
    
    if (plot_counts && plot_norm_expression) {
      counts_for_plot <- tryCatch({
        reshape2::melt(counts_matrix_named[counts_matrix_named$hgnc_symbol == genes,])
      }, error = function(e) {
        stop('Error in accessing data. Counts matrix must be stored as counts_matrix_named to plot counts with this function. Gene symbols must be a column called hgnc_symbol.')
      })
      norm_counts_for_plot <- tryCatch({
        reshape2::melt(vst_norm_matrix_named[vst_norm_matrix_named$hgnc_symbol == genes,])
      }, error = function(e) {
        stop('Error in accessing data. Normalized expression matrix must be stored as vst_norm_matrix_named to plot normalized expression with this function. Gene symbols must be a column called hgnc_symbol.')
      })
      
      if (nrow(counts_for_plot) == 0) {
        stop(paste(genes, 'not present in dataset.'))
      }
      
      counts_for_plot$group <- substr(counts_for_plot$variable, 1,2)
      
      # plot
      counts_plot <- ggplot(counts_for_plot %>% filter(hgnc_symbol == genes), aes(x = group, y = value, fill = group)) +
        geom_boxplot(size = 0.75) +
        scale_fill_manual(values = c('#CCCCCC', '#EF778E')) +
        geom_point(data = norm_counts_for_plot %>% filter(hgnc_symbol == genes), aes(x = group, y = value), size = 2.5) +
        theme_classic() +
        ylab('counts') +
        ggtitle(gene, 'raw counts') +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = 20, color = "black"),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.line = element_line(size = 0.75),
              axis.ticks = element_line(size = 0.75),
              plot.title = element_text(size = 18, hjust = 0.5))
      
      norm_plot <- ggplot(norm_counts_for_plot %>% filter(hgnc_symbol == genes), aes(x = group, y = value, fill = group)) +
        geom_boxplot(size = 0.75) +
        scale_fill_manual(values = c('#CCCCCC', '#EF778E')) +
        geom_point(data = norm_counts_for_plot %>% filter(hgnc_symbol == genes), aes(x = group, y = value), size = 2.5) +
        theme_classic() +
        ylab('normalized counts') +
        ggtitle(gene, 'normalized expression') +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = 20, color = "black"),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.line = element_line(size = 0.75),
              axis.ticks = element_line(size = 0.75),
              plot.title = element_text(size = 18, hjust = 0.5))
      
      plots <- plot_grid(counts_plot, norm_plot, nrow = 1, align = "hv", axis = "lrbt")
      return(plots)
      
    } else if (!plot_counts && plot_norm_expression) {
      if (plot_counts) {
        counts_for_plot <- tryCatch({
          reshape2::melt(counts_matrix_named[counts_matrix_named$hgnc_symbol == genes,])
        }, error = function(e) {
          stop('Error in accessing data. Counts matrix must be stored as counts_matrix_named to plot counts with this function. Gene symbols must be a column called hgnc_symbol.')
        })
      } else if (plot_norm_expression) {
        counts_for_plot <- tryCatch({
          reshape2::melt(vst_norm_matrix_named[vst_norm_matrix_named$hgnc_symbol == genes,])
        }, error = function(e) {
          stop('Error in accessing data. Normalized expression matrix must be stored as vst_norm_matrix_named to plot normalized expression with this function. Gene symbols must be a column called hgnc_symbol.')
        }) 
      } else {
        stop('One of plot_counts or plot_norm_expression must be TRUE.')
      }
      
      if (nrow(counts_for_plot) == 0) {
        stop(paste(genes, 'not present in dataset.'))
      }
      
      counts_for_plot$group <- substr(counts_for_plot$variable, 1,2)
      
      # plot
      pval <- dp_res_df_named[dp_res_df_named$hgnc_symbol == genes,]$padj %>% round(., digits = 4)
      plot <- ggplot(counts_for_plot, aes(x = group, y = value, fill = group)) +
        geom_boxplot(size = 0.75) +
        scale_fill_manual(values = c('#CCCCCC', '#EF778E')) +
        geom_point(data = counts_for_plot, aes(x = group, y = value), size = 2.5) +
        theme_classic() +
        ylab('normalized counts') +
        ggtitle(paste(genes, 'P =', pval)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = 20, color = "black"),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.line = element_line(size = 0.75),
              axis.ticks = element_line(size = 0.75),
              plot.title = element_text(size = 18, hjust = 0.5))
      
      
    }  
  } else {
    
    ## For multiple genes - only plots normalized counts for now (will update in future)
    
    norm_counts_for_plot <- reshape2::melt(vst_norm_matrix_named[vst_norm_matrix_named$hgnc_symbol %in% genes,])
    norm_counts_for_plot$group <- substr(norm_counts_for_plot$variable, 1, 2)
    
    plotlist <- NULL
    for (gene in genes) {
      pval <- dp_res_df_named[dp_res_df_named$hgnc_symbol == gene,]$padj %>% round(., digits = 4)
      plotlist[[gene]] <- ggplot(norm_counts_for_plot %>% filter(hgnc_symbol == gene), aes(x = group, y = value, fill = group)) +
        geom_boxplot(size = 0.75) +
        scale_fill_manual(values = c('#CCCCCC', '#EF778E')) +
        geom_point(data = norm_counts_for_plot %>% filter(hgnc_symbol == gene), aes(x = group, y = value), size = 2.5) +
        theme_classic() +
        ylab('normalized counts') +
        ggtitle(paste(gene, 'P =', pval)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = 20, color = "black"),
              axis.text.y = element_text(size = 18, color = "black"),
              axis.line = element_line(size = 0.75),
              axis.ticks = element_line(size = 0.75),
              plot.title = element_text(size = 18, hjust = 0.5))
    }
    plot <- plot_grid(plotlist = plotlist, align = "hv", axis = "lrbt")
    
  }
  return(plot)
}