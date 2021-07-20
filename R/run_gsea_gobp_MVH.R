#' @title run_gsea_gobp_MVH
#' @description Run fGSEA on GO biological processes pathway set.
#' @param marker_list FindMarkers output, not converted to human.
#' @param filter_pval if T, filters FindMarkers results to remove p > 0.05.
#' @export
#' @return fGSEA results as dataframe.

run_gsea_gobp_MVH <- function(marker_list, filter_pval = T) {
  
  require(msigdbr)
  require(fgsea)
  require(mvhfunctions)
  
  # check for msigdbr loaded already
  message('Loading pathways...\n')
  human.genes <- msigdbr(species = "Homo sapiens")
  genesets.GOBP <- filter(human.genes, gs_subcat == "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name)
  
  # filter for significant p vals if selected above
  if (filter_pval) {
    marker_list <- marker_list %>% filter(p_val_adj <= 0.05)
  }
  
  message('Converting genes to human...\n')
  # convert to human
  marker_list_human <- marker_list %>% convert_FindMarkers_human_MVH()
  
  message('Running fGSEA...\n')
  
  # run gsea
  set.seed(3)
  gsea_gobp <- fgsea(pathways = genesets.GOBP,
                     stats = marker_list_human,
                     nproc = 1) %>%
    arrange(desc(NES)) %>% 
    mutate(log10pval = -log10(pval)) %>%
    select(pathway, NES, pval, log10pval, everything())
  
  message('Done!\n')
  return(gsea_gobp)
  
}