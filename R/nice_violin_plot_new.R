#' @title nice_violin_plot
#' @description nicer looking version of the Seurat function VlnPlot.
#' @param seurat_obj Seurat object.
#' @param features genes to plot.
#' @param group_by how to group your data, i.e the X axis of the plot - usually a metadata column in the seurat object.
#' @param cols optional: colours to fill the violins with.
#' @param pt.size size of plotted points.
#' @param sort sort X axis in descending order of average expression?
#' @param n_col if multiple features, how many columns to create in the final plot_grid object.
#' @param plot_hline plot a dashed horizontal line at 0?
#' @param perform_stats perform pairwise wilcoxon rank sum test between groups and print P values? *** only works for single plots for now
#' @export
#' @return violin plot.
#'

nice_violin_plot <- function(seurat_obj, features, group_by = NULL, cols = NULL, pt.size = 0.3, sort = T, n_col = NULL, plot_hline = T, perform_stats = T) {
  
  require(Seurat)
  require(tidyverse)
  require(pals)
  require(reshape2)
  
  # determine number of violins
  if (!is.null(group_by)) {
    groups <- seurat_obj[[]] %>% dplyr::select(all_of(group_by)) %>% unique()
    n_groups <- nrow(groups)
  } else {
    groups <- as.character(unique(Idents(seurat_obj)))
    n_groups <- length(groups)
  }
  
  # set colormap 
  if ((is.null(cols)) & (n_groups < 13)) {
    cols <- tol(n_groups)
  }
  
  if (length(features) == 1) {
    plot <- Seurat::VlnPlot(seurat_obj,
                            group.by = group_by,
                            features = features,
                            pt.size = pt.size,
                            cols = cols,
                            sort = sort) +
      Seurat::NoLegend() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))
    
    if (perform_stats) {
      
      if (pt.size == 0) {
        message('pt.size cannot be 0 if performing stats.')
      } else {
        ## perform stats: only works on single plots for now ##
        plot_data_all <- ggplot_build(plot)
        x_groups <- plot_data_all$layout$panel_params[[1]]$x$get_labels() # get X groups
        plot_data <- plot_data_all$data[[2]] # expression data per cell
        
        # add information about groups to plot_data
        plot_data$group_name <- NA
        for (ii in 1:n_groups) {
          plot_data[plot_data$group == ii,]$group_name <- x_groups[ii]
        }
        
        # perform stats on plot_data
        if (n_groups > 2) {
          stats <- pairwise.wilcox.test(plot_data$y, plot_data$group_name, p.adjust.method = "bonferroni")
          p_vals <- stats$p.value %>% data.frame() %>% rownames_to_column(var = "group1") %>% melt()
          p_vals <- p_vals[!is.na(p_vals$value),]
          p_vals$value <- signif(p_vals$value, digits = 3)
          
          for (ii in 1:nrow(p_vals)) {
            message(paste(p_vals$group1[ii], 'vs', p_vals$variable[ii], 'P =', p_vals$value[ii]))
          }
          
        } else if (n_groups == 2) {
          stats <- wilcox.test(plot_data %>% filter(group_name == x_groups[1]) %>% pull(y),
                               plot_data %>% filter(group_name == x_groups[2]) %>% pull(y))
          pval <- signif(stats$p.value, digits = 3)
          message(paste(x_groups[1], 'vs', x_groups[2], 'P =', pval))
          # add pval to plot
          plot <- plot + 
            labs(title = features, subtitle = paste('P = ', pval)) +
            theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5))
        }
      }
    }
    
    if (plot_hline == T) {
      plot <- plot + geom_hline(yintercept = 0, linetype = "dashed")
      return(plot)
    } else {
      return(plot)
    }
  } else {
    
    if (perform_stats) {
      message('Stats cannot be performed if plotting multiple features (for now).')
    }
    
    plotlist <- Seurat::VlnPlot(seurat_obj,
                                group.by = group_by,
                                features = features,
                                pt.size = pt.size,
                                sort = sort,
                                cols = cols,
                                combine = F)
    plotlist <- lapply(plotlist, function(x)
      x + theme(axis.title.x = element_blank(),
                axis.title.y = element_text(size = 10)) + NoLegend())
    
    if (plot_hline == T) {
      plotlist <- lapply(plotlist, function(x)
        x + geom_hline(yintercept = 0, linetype = "dashed"))
    }
    plots <- cowplot::plot_grid(plotlist = plotlist, ncol = n_col)
    return(plots)
  }
}
