#' @title add_gene_names_deseq2
#' @description Add gene names to deseq2 results, counts matrix or normalized expression matrix. Can only do one at a time for now. 
#' @param deseq_results deseq2 results data frame.
#' @param deseq_norm_matrix normalized expression matrix.
#' @param counts_matrix counts matrix. 
#' @export
#' @return named expression data.



add_gene_names_deseq2 <- function(deseq_results = NULL, deseq_norm_matrix = NULL, counts_matrix = NULL) {
  
  require(biomaRt)
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  
  if (!is.null(deseq_results)) {
    genes <- deseq_results$Ensembl
  } else if (!is.null(deseq_norm_matrix)) {
    genes <- rownames(deseq_norm_matrix)
  } else if (!is.null(counts_matrix)) {
    genes <- counts_matrix$Ensembl
  }
  
  data <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = genes,
                mart = ensembl)
  
  # DESeq2 results
  if (!is.null(deseq_results)) {
    deseq_results_named <- merge(x = deseq_results,
                                 y = data,
                                 by.x = "Ensembl",
                                 by.y = "ensembl_gene_id") %>% dplyr::select(Ensembl, hgnc_symbol, everything())
    
    # for genes without gene name, add ENSEMBL ID to the gene name column.
    deseq_results_named$hgnc_symbol <- case_when(
      deseq_results_named$hgnc_symbol == "" ~ deseq_results_named$Ensembl,
      TRUE ~ deseq_results_named$hgnc_symbol
    )
    return(deseq_results_named)
  }
 
  # counts matrix
  if (!is.null(counts_matrix)) {
    counts_matrix <- counts_matrix %>% rownames_to_column(var = "Ensembl")
    counts_matrix_named <- merge(x = counts_matrix,
                                y = data,
                                by.x = "Ensembl",
                                by.y = "ensembl_gene_id") %>% dplyr::select(Ensembl, hgnc_symbol, everything())
    
    # for genes without gene name, add ENSEMBL ID to the gene name column.
    counts_matrix_named$hgnc_symbol <- case_when(
      counts_matrix_named$hgnc_symbol == "" ~ counts_matrix_named$Ensembl,
      TRUE ~ counts_matrix_named$hgnc_symbol
    )
    return(counts_matrix_named)
  }
  
  # normalized matrix
  if (!is.null(deseq_norm_matrix)) {
    deseq_norm_matrix <- data.frame(deseq_norm_matrix) %>% rownames_to_column(var = "Ensembl")
    deseq_norm_matrix_named <- merge(x = deseq_norm_matrix,
                                   y = data,
                                   by.x = "Ensembl",
                                   by.y = "ensembl_gene_id") %>% dplyr::select(Ensembl, hgnc_symbol, everything())
   
    # for genes without gene name, add ENSEMBL ID to the gene name column.
    deseq_norm_matrix_named$hgnc_symbol <- case_when(
      deseq_norm_matrix_named$hgnc_symbol == "" ~ deseq_norm_matrix_named$Ensembl,
      TRUE ~ deseq_norm_matrix_named$hgnc_symbol
    )
    return(deseq_norm_matrix_named)
  }
}