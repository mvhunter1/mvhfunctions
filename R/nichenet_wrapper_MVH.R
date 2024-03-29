#' @title nichenet_wrapper_MVH
#' @description perform Nichenet analysis on clusters of interest.
#' @param seurat_obj Seurat object.
#' @param sender_ident Name of cluster you want to be considered as the "sender".
#' @param receiver_ident Name of cluster you want to be considered as the "receiver".
#' @param ligand_target_cutoff Cutoff for ligand-target analysis.
#' @export
#' @return Nichenet results plots.

nichenet_wrapper_MVH <- function(seurat_obj, sender_ident, receiver_ident, ligand_target_cutoff = 0.5) {
  
  require(Seurat)
  require(tidyverse)
  require(cowplot)
  require(nichenetr)
  
  # load necessary files
  setwd('/Volumes/GoogleDrive-107501420737632855873/.shortcut-targets-by-id/1L5uwxnNO44Pd7GBMBhMGO9vvsi7uifYZ/White Lab Gsuite Main Drive/Lab members/Miranda Hunter/Miranda_R_new/nichenet/')
  load('nichenet_ligand_target_matrix.R')
  load('nichenet_LR_network.R')
  load('nichenet_weighted_networks.R')
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  
  message("Determining expressed genes in sender and receiver cells...")
  expressed_genes_receiver <- get_expressed_genes(ident = receiver_ident, 
                                                  seurat_obj = seurat_obj,
                                                  assay_oi = "SCT",
                                                  pct = 0.1)
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  expressed_genes_sender <- get_expressed_genes(ident = sender_ident, 
                                                seurat_obj = seurat_obj, 
                                                assay_oi = "SCT",
                                                pct = 0.1) %>% unlist() %>% unique()
  
  message("Calculating upregulated genes in receiver cells...")
  DE_table_receiver <- FindMarkers(object = seurat_obj, 
                                   ident.1 = receiver_ident, 
                                   min.pct = 0.1,
                                   assay = "SCT") %>% rownames_to_column(var = "gene")
  
  geneset_oi <- DE_table_receiver %>% 
    filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% 
    pull(gene)
  
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  message("Determining which ligands and receptors are expressed in the sender and receiver cells...")
  ligands <- lr_network %>% pull(from) %>% unique()
  receptors <- lr_network %>% pull(to) %>% unique()
  
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  
  message("Determining potential ligands and predicting ligand activities...")
  potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  ligand_activities <- predict_ligand_activities(geneset_oi,
                                                 background_expressed_genes,
                                                 ligand_target_matrix,
                                                 potential_ligands)
  ligand_activities <- ligand_activities %>% 
    arrange(-pearson) %>% 
    mutate(rank = rank(desc(pearson)))
  
  best_upstream_ligands <- ligand_activities %>% 
    top_n(20, pearson) %>% 
    arrange(-pearson) %>% 
    pull(test_ligand)
  
  message("Calculating ligand targets...")
  active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                    geneset = geneset_oi,
                                                                    ligand_target_matrix,
                                                                    n = 250) %>% bind_rows() %>% drop_na()
  
  message("Preparing data for plotting...")
  active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                    ligand_target_matrix,
                                                                    cutoff = ligand_target_cutoff) # changed to top 10% for integrated object
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  
  order_targets <- active_ligand_target_links_df$target %>% 
    unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% 
    make.names()
  
  vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()
  
  ligand_pearson_matrix <- ligand_activities %>% 
    dplyr::select(pearson) %>% 
    as.matrix() %>% 
    magrittr::set_rownames(ligand_activities$test_ligand)
  
  vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>% 
    as.matrix(ncol = 1) %>% 
    magrittr::set_colnames("PCC")
  
  lr_network_top <- lr_network %>% 
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  
  message("Determining potential receptors expressed on receiver cells...")
  best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()
  lr_network_top_df_large <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df <- lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors <- dist(lr_network_top_matrix, method = "binary")
  hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
  order_receptors <- hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()
  
  message("Generating plots...")
  # plots
  # ligand activity predictions
  ligand_activity_preds <- vis_ligand_pearson %>% make_heatmap_ggplot(paste0("ligands (", sender_ident, ")"), "predicted ligand activity", 
                                                                      color = "darkorange", 
                                                                      legend_position = "right", 
                                                                      x_axis_position = "top", 
                                                                      legend_title = "Pearson\ncorrelation\ncoefficient") + 
    theme(axis.title = element_text(face = "bold"), 
          axis.text = element_text(size = 12, angle = 0), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 10), 
          axis.text.x.top = element_blank(),
          axis.title.x.top = element_text(size = 14), 
          legend.title.align = 0.5, 
          axis.title.y.left = element_text(size = 14, margin = margin(r = 20, l = 15)), 
          axis.text.y.left = element_text(margin = margin(r = 0))) 
  
  
  # predicted receptors
  receptor_preds <- vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot(paste0("ligands (", sender_ident, ")"), paste0("predicted receptors (", receiver_ident, ")"), 
                                                                                color = "mediumvioletred", 
                                                                                x_axis_position = "top", 
                                                                                legend_title = "interaction\npotential") + 
    theme(axis.title = element_text(size = 14, face = "bold"), 
          axis.text.y.left = element_text(size = 12), 
          legend.title = element_text(size = 10, face = "bold", hjust = 0.5), 
          axis.text.x.top = element_text(angle = 90, vjust = 0.5, size = 12), 
          axis.title.x.top = element_text(margin = margin(b = 12)), 
          axis.title.y.left = element_text(margin = margin(r = 10, l = 60)),
          legend.position = "right")
  
  first_row <- plot_grid(ligand_activity_preds, receptor_preds, rel_widths = c(1,3)) 
  
  # predicted targets
  target_preds <- vis_ligand_target %>% make_heatmap_ggplot(paste0("ligands (", sender_ident, ")"), paste0("predicted target genes (", receiver_ident, ")"),
                                                            color = "purple",
                                                            x_axis_position = "top",
                                                            legend_title = "regulatory\npotential") + 
    scale_fill_gradient2(low = "whitesmoke", high = "purple", breaks = c(0,0.005,0.01)) + 
    theme(axis.title = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"), 
          axis.title.y = element_text(margin = margin(l = 30, r = 10)), 
          axis.title.x.top = element_text(margin = margin(b = 10, t = 50)), 
          axis.text.x.top = element_text(vjust = 0.5, size = 4), 
          legend.position = "right", 
          legend.margin = margin(l = 10, r = 25),
          legend.text = element_text(size = 10),
          axis.text.y.left = element_text(size = 12),
          plot.margin = margin(r = 20, b = 20))
  
  
  plots <- plot_grid(first_row, target_preds, nrow = 2)
  message("Done!")
  return(plots)
}