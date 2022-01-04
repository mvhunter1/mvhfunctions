#' @title invasion_assay_plotting_function
#' @description Plot normalized data from invasion assay.
#' @param plot_data data frame with normalized results.
#' @param sep_by_replicates Plot each replicate separately?
#' @export
#' @return Plots with results from the experiment. If sep_by_replicates = F and there are multiple replicates, they will be labelled with different shapes on the plot.



invasion_assay_plotting_function <- function(plot_data, sep_by_replicates = F) {
  
  
  if ("replicate" %in% colnames(plot_data)) {
    plot_data$replicate <- as.character(plot_data$replicate)
    replicates <- sort(unique(plot_data$replicate))
    n_replicates <- length(replicates)
    plot_data$replicate <- factor(plot_data$replicate, levels = replicates)
    
    if (sep_by_replicates) {
      
      plotlist <- NULL
      for (n_rep in replicates) {
        rep_plot_data <- plot_data %>% filter(plot_data$replicate == n_rep) 
        plotlist[[n_rep]] <- ggplot(rep_plot_data, aes(x = group, y = num_cells)) +
          geom_boxplot() +
          geom_point(size = 3) +
          geom_hline(yintercept = 1, linetype = "dashed") +
          theme(axis.title.x = element_blank(),
                axis.text = element_text(size = 18, color = "black"),
                axis.title = element_text(size = 20, color = "black")) +
          ylab("relative # of invaded cells")
        
      }
      plot <- plot_grid(plotlist = plotlist, ncol = n_replicates)
      
    } else {
      plot <- ggplot(plot_data, aes(x = group, y = num_cells)) +
        geom_boxplot() +
        geom_jitter(aes(x = group, y = num_cells, color = replicate), size = 3, width = 0.1) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        theme(axis.title.x = element_blank(),
              axis.text = element_text(size = 18, color = "black"),
              axis.title = element_text(size = 20, color = "black")) +
        ylab("relative # of invaded cells")
    }
    
    
  } else {
    plot <- ggplot(plot_data, aes(x = group, y = num_cells)) +
      geom_boxplot() +
      geom_point(size = 3) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      theme(axis.title.x = element_blank(),
            axis.text = element_text(size = 18, color = "black"),
            axis.title = element_text(size = 20, color = "black")) +
      ylab("relative # of invaded cells")
  }
  return(plot)
}