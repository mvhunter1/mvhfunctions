#' @title process_cyquant_plate
#' @description Process and plot data from a 96 well plate used with Cyquant assay. This is only set up to use with a specific assay and plate layout - not for general use.
#' @param path_to_cyquant_plate Path to the excel file with the data from the plate reader.
#' @export
#' @return Plots with results from the experiment.


process_cyquant_plate <- function(path_to_cyquant_plate) {
  
  require(tidyverse)
  require(readxl)
  
  cyquant_plate <- readxl::read_excel(path_to_cyquant_plate) %>%
    column_to_rownames(var = "...1") %>%
    dplyr::select(-...14)
  
  # remove the last 2 rows since they don't contain data
  cyquant_plate <- cyquant_plate[1:6,]
  
  cell_conc <- c("2500", "5000", "10,000", "20,000")
  
  cyquant_plate_norm <- NULL
  
  for (cell_num in cell_conc) {
    
    index <- which(cell_conc == cell_num)
    first_col <- 1 + 3*(index - 1) # find columns of matrix with relevant data
    last_col <- first_col + 2
    
    data <- cyquant_plate[,first_col:last_col] %>% as.matrix()
    mean_control <- mean(data[1,])
    data_norm <- data / mean_control
    
    if (is.null(cyquant_plate_norm)) {
      cyquant_plate_norm <- data_norm
    } else {
      cyquant_plate_norm <- cbind(cyquant_plate_norm, data_norm)
    }
  }
  
  hmgb2_conc <- c("0", "5", "10", "30", "100", "250")
  
  plots <- NULL
  
  for (cell_num in cell_conc) {
    
    index <- which(cell_conc == cell_num)
    first_col <- 1 + 3*(index - 1) # find columns of matrix with relevant data
    last_col <- first_col + 2
    
    data <- cyquant_plate_norm[,first_col:last_col] 
    
    # reorganize data for plotting
    data_plot_organized <- NULL
    for (conc in hmgb2_conc) {
      index_conc <- which(hmgb2_conc == conc)
      data_plot <- data.frame(conc_hmg = rep(conc, 3),
                              num_cells = data[index_conc,])
      
      if (is.null(data_plot_organized)) {
        data_plot_organized <- data_plot
      } else {
        data_plot_organized <- rbind(data_plot_organized, data_plot)
      }
    }
    
    data_plot_organized$conc_hmg <- factor(data_plot_organized$conc_hmg, levels = hmgb2_conc)
    plots[[index]] <- ggplot(data_plot_organized, aes(x = conc_hmg, y = num_cells)) +
      geom_boxplot(width = 0.65) +
      geom_point(size = 2.5) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      xlab("ng/mL purified HMGB2") +
      ylab("normalized cell number") +
      ggtitle(paste(cell_num, "cells/well")) +
      #coord_cartesian(ylim = c(0.5, 4)) +
      theme(axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 18, color = "black"),
            plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
            axis.title.y = element_text(margin = margin(l = 20))) 
    
  }
  
  plots <- plot_grid(plotlist = plots, ncol = 4)
  return(plots)
}