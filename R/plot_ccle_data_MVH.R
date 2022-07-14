#' @title plot_ccle_data_MVH
#' @description Plot RNA-seq data from CCLE/DepMap repository. ccle_data.R, aka CCLE_expression.csv, and ccle_sample_info.R, aka sample_info.csv, must be loaded first for this function to work. Both files can be downloaded from: https://depmap.org/portal/download/
#' @param ccle_data CCLE log2 RSEM data, with genes as rownames as HGNC ID.
#' @param genes genelist to plot.
#' @param cell_line_id Cell line you want to plot. If A375 or SKMEL28, can just put those names. If anything else, need to put DepMap cell line ID.
#' @export
#' @return bargraph of normalized expression of the given genes.


plot_ccle_data_MVH<- function(genes, cell_line_id = "A375") {
  # if not A375 or SKMEL28, need to provide cell_line_id as DepMap ID.
  
  load('/Volumes/GoogleDrive-107501420737632855873/.shortcut-targets-by-id/1L5uwxnNO44Pd7GBMBhMGO9vvsi7uifYZ/White Lab Gsuite Main Drive/Lab members/Miranda Hunter/Miranda_R_new/CCLE/ccle_data.R')
  load('/Volumes/GoogleDrive-107501420737632855873/.shortcut-targets-by-id/1L5uwxnNO44Pd7GBMBhMGO9vvsi7uifYZ/White Lab Gsuite Main Drive/Lab members/Miranda Hunter/Miranda_R_new/CCLE/ccle_sample_info.R')
  
  if (cell_line_id == "A375") {
    cell_line_id <- "ACH-000219"
  } else if (cell_line_id == "SKMEL28") {
    cell_line_id <- "ACH-000615"
  } else {
    # check if cell line ID is in dataset. if not, throw an error.
    if (!cell_line_id %in% rownames(ccle_data)) {
      stop('Cell line ID is not found in dataset.')
    }
  }
  
  cell_line_data <- ccle_data[rownames(ccle_data) %in% cell_line_id,] %>% 
    t() %>% 
    data.frame() %>% 
    rownames_to_column(var = "gene")
  colnames(cell_line_data)[2] <- "gene_exp"
  
  # Check if any genes of interest are not expressed in CCLE dataset.
  genes_not_expressed <- setdiff(genes, cell_line_data$gene)
  if (length(genes_not_expressed) > 0) {
    message(paste(genes_not_expressed, "is not expressed in CCLE data."))
  }
  
  # Plot.
  plot_data <- cell_line_data[cell_line_data$gene %in% genes,] %>% arrange(-gene_exp)
  plot_data$gene <- factor(plot_data$gene, levels = plot_data$gene)
  
  plot <- ggplot(plot_data, aes(x = gene, y = gene_exp, fill = gene_exp)) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = rev(viridis(n = 100))) +
    theme(axis.text = element_text(size = 16, color = "black"),
          axis.text.x = element_text(angle = -45, hjust = 0.2, vjust = 0),
          axis.title.y = element_text(size = 18, color = "black", margin = margin(r = 10)),
          axis.title.x = element_blank(),
          legend.position = "none") +
    ylab("normalized expression\n(log2 RSEM)")
  
  # Use sample_info.csv to get cell line names.
  cell_line_name <- ccle_sample_info[ccle_sample_info$DepMap_ID == cell_line_id,]$stripped_cell_line_name
  # Add cell line name to plot.
  plot <- plot + ggtitle(paste("Expression in", cell_line_name)) + theme(plot.title = element_text(size = 24, hjust = 0.5))
  return(plot)
    
}