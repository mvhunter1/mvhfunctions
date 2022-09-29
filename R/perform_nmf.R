#' @title perform_nmf
#' @description Perform NMF and create gene modules from expression matrix. Slightly modified code from Reuben and Dalia.
#' @param seurat_obj Seurat object.
#' @param rank Rank for NMF, default = 20. 
#' @param cluster_name Name of cluster or clusters you want to run NMF on.
#' @export
#' @return NMF results list with matrix, modules, scores and coef.

perform_nmf <- function(seurat_obj, rank = 20, cluster_name) {
  
  require(Seurat)
  require(RcppML) 
  require(Matrix)
  
  seurat_cluster <- subset(seurat_obj, idents = cluster_name)
  seurat_cluster_mat <- as.matrix(GetAssayData(seurat_cluster))
  
  seurat_cluster <- FindVariableFeatures(seurat_cluster, selection.method = "vst", verbose = F)
  seurat_cluster_var_genes <- VariableFeatures(seurat_cluster)
  
  expression_matrix <- seurat_cluster_mat[rownames(seurat_cluster_mat) %in% seurat_cluster_var_genes,] 
  
  # make sure matrix is stored as matrix
  expression_matrix <- as.matrix(expression_matrix)
  
  # check to see if there are any negative values and if so, set to 0
  expression_matrix[expression_matrix < 0] <- 0
  
  # remove any rows where expression = 0, otherwise it will give an error when running the NMF function
  row_sums <- rowSums(expression_matrix)
  zero_rows <- which(row_sums > 0)
  expression_matrix_filt <- expression_matrix[zero_rows,]
  
  # run NMF
  message('Running NMF...')
  t1 <- Sys.time() 
  
  if (length(rank) > 1) {
    stop('Multiple ranks cannot be used with use_fast_method.')
  }
  
  nmf_results <- RcppML::nmf(A = Matrix::Matrix(expression_matrix_filt, sparse = T), 
                             k = rank, 
                             seed = 'ica') 
  
  t2 <- Sys.time() 
  message('NMF calculations complete!')
  difftime(t2,t1)
  
  # Make modules
  message('Making modules...')
  
  scores <- nmf_results$w
  coef <- nmf_results$h
  rownames(scores) <- rownames(expression_matrix_filt) # critical!!!!
  colnames(coef)<- colnames(expression_matrix_filt) # critical!!!!
  
  
  o <- order(apply(scores, 2, var), decreasing = TRUE)
  scores <- scores[, o]
  coef <- coef[o, ]
  
  # Make modules
  scores.norm <- t(scores) / apply(scores, 2, mean)
  ranks_x <- t(apply(-t(scores.norm), 1, base::rank))
  ranks_y <- apply(-t(scores.norm), 2, base::rank)
  
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  
  modules <- apply(ranks_y, 2, function(m){
    a <- sort(m[is.finite(m)])
    a <- a[a == 1:length(a)]
    names(a)
  })
  
  l <- sapply(modules, length)
  keep <- (l >= 5) # remove any module with less than 5 genes
  scores <- scores[, keep]
  coef <- coef[keep, ]
  
  # Repeat module making
  ranks_x <- t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y <- apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules <- apply(ranks_y, 2, function(m){
    a <- sort(m[is.finite(m)])
    a <- a[a == 1:length(a)]
    names(a)
  })
  
  # store results in list
  results_list <- list(results = nmf_results, coef = coef, scores = scores, modules = modules)
  
  message('Done!')
  return(results_list)
 
}