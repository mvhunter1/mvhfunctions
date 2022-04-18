#' @title go_enrichment_nmf
#' @description Perform GO term enrichment on NMF results using GProfiler. Both GO BP and GO CC enrichment will be run.
#' @param nmf_results nmf_results, output from perform_nmf.
#' @param num_genes number of top scoring genes per factor to be used for GO analysis.
#' @param n_terms_plot number of GO terms to be plotted per factor.
#' @param species species to use for GO enrichment. must be one of 'fish', 'zebrafish', or 'human'.
#' @export
#' @return results as list including dataframes with all results, and plots. 

go_enrichment_nmf <- function(nmf_results, num_genes = 150, n_terms_plot = 10, species = 'fish') {
  
  require(gprofiler2)
  set.seed(100)
  species <- tolower(species)
  
  if (species %in% c('fish', 'zebrafish')) {
    message('Zebrafish genes will be converted to human before GO analysis.')
  } else if (species == 'human') {
    message('Human genesets will be used for GO analysis.')
  } else if (!species %in% c('fish', 'zebrafish', 'human')) {
    stop('Species must be one of fish, zebrafish, or human.')
  }
  
  scores <- nmf_results$scores
  BP_plots <- NULL
  CC_plots <- NULL
  BP_results_allmodules <- NULL
  CC_results_allmodules <- NULL
  
  require(mvhfunctions)
  
  for (ii in 1:dim(scores)[2]){
    
    message(paste0('Performing GO enrichment on factor ', ii, '...'))
    # Retrieve top scoring genes per NMF factor
    query_genes <- rownames(scores)[order(scores[,ii], decreasing = TRUE)[1:num_genes]] 
    
    if (species %in% c('fish', 'zebrafish')) {
      require(mvhfunctions)
      query_genes <- convert_to_human_list(query_genes)
    }
    
    # Perform GO enrichment
    enrichment <- gost(query_genes,
                       organism = 'hsapiens',
                       sources = 'GO',
                       significant = T)
    
    # If no enrichment, continue loop to next iteration
    if (is.null(enrichment)){
      message(paste0('No significant pathways found for NMF factor '), ii, '...')
      next
    }
    
    # summarize results for plotting
    enrichment_res_BP <- enrichment$result %>% 
      filter(source == "GO:BP") %>%
      mutate(log10pval = -log10(p_value)) %>%
      arrange(-log10pval) %>% mutate(NMF_module = ii)
    
    enrichment_res_CC <- enrichment$result %>% 
      filter(source == "GO:CC") %>%
      mutate(log10pval = -log10(p_value)) %>%
      arrange(-log10pval) %>% mutate(NMF_module = ii)
    
    # organize data for plotting
    enrichment_res_BP$term_name <- factor(enrichment_res_BP$term_name, levels = rev(enrichment_res_BP$term_name))
    enrichment_res_CC$term_name <- factor(enrichment_res_CC$term_name, levels = rev(enrichment_res_CC$term_name))
    
    n_terms_BP <- ifelse(n_terms_plot > nrow(enrichment_res_BP), nrow(enrichment_res_BP), n_terms_plot)
    n_terms_CC <- ifelse(n_terms_plot > nrow(enrichment_res_CC), nrow(enrichment_res_CC), n_terms_plot)
    
    BP_plots[[ii]] <- ggplot(enrichment_res_BP[1:n_terms_BP,], aes(x = log10pval, y = term_name)) +
      geom_bar(stat = "identity", fill = '#8ECED5') +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 13, color = "black"),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14, color = "black"),
            axis.title.x = element_text(size = 16, color = "black"),
            plot.title = element_text(size = 18, color = "black", hjust = 0.5)) +
      xlab('-log10 p-value') +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
      ggtitle(paste('GO biological processes: NMF factor', ii))
    
    
    CC_plots[[ii]] <- ggplot(enrichment_res_CC[1:n_terms_CC,], aes(x = log10pval, y = term_name)) +
      geom_bar(stat = "identity", fill = '#CED58E') +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 13, color = "black"),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14, color = "black"),
            axis.title.x = element_text(size = 16, color = "black"),
            plot.title = element_text(size = 18, color = "black", hjust = 0.5)) +
      xlab('-log10 p-value') +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
      ggtitle(paste('GO cellular component: NMF factor', ii))
    
    # # plot expression of top genes per factor
    # genes <- rownames(scores)[order(scores[,ii],decreasing = TRUE)[1:num_genes]]
    # exp <- colMeans(GetAssayData(EF.filt, slot = 'scale.data')[genes,colnames(fish_tumor)])
    # EF.filt <- AddMetaData(EF.filt, metadata = exp, col.name = paste0('NMF_', ii, '_expr'))
    
    if (is.null(BP_results_allmodules)) {
      BP_results_allmodules <- enrichment_res_BP
    } else {
      BP_results_allmodules <- rbind(BP_results_allmodules, enrichment_res_BP)
    }
    
    if (is.null(CC_results_allmodules)) {
      CC_results_allmodules <- enrichment_res_CC
    } else {
      CC_results_allmodules <- rbind(CC_results_allmodules, enrichment_res_CC)
    }
    
  }
  
  results <- list(BP_results_all = BP_results_allmodules,
                  CC_results_all = CC_results_allmodules,
                  BP_plots = BP_plots,
                  CC_plots = CC_plots)
  return(results)
}