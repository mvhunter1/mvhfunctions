#' @title perform_nmf
#' @description perform NMF on a matrix of expression data.
#' @param expression_matrix For scRNA-seq or Visium, typically the output from GetAssayData for your cluster/cell type of interest.
#' @param rank Rank used for NMF, can be determined using NMFestimaterank.
#' @param save_NMF_results if T, will save NMF function output as .R file. This can be useful since the function takes several hours to run.
#' @param filename if save_NMF_results = T, filename that the NMF results will be saved as.
#' @export
#' @return a list containing the NMF results matrix, coefficients, scores, and modules.

perform_nmf <- function(expression_matrix, rank = 20, save_NMF_results = F, filename = NULL) {
  
  # modified code from Reuben
  
  require(NMF)
  require(fastICA)
  
  # make sure matrix is saved as matrix
  expression_matrix <- as.matrix(expression_matrix)
  
  # check to see if there are any negative values and if so, set to 0
  expression_matrix[expression_matrix < 0] <- 0
  
  # remove any rows where expression = 0, otherwise it will give an error when running the NMF function
  row_sums <- rowSums(expression_matrix)
  zero_rows <- which(row_sums > 0)
  expression_matrix_filt <- expression_matrix[zero_rows,]
  
  # run NMF
  message('Running NMF...')
  nmf_results <- nmf(x = expression_matrix_filt, rank = rank, seed = 'ica', method = 'nsNMF')
  message('NMF calculations complete!')
  
  if (save_NMF_results) {
    message('Saving NMF results...')
    if (is.null(filename)) {
      save(nmf_results, file = paste0("nmf_results_", Sys.time(), ".R"))
    } else {
      save(nmf_results, file = filename)
    }
  }
  
  # Make modules
  message('Making modules...')
  scores <- basis(nmf_results)
  coef <- coefficients(nmf_results)
  
  o <- order(apply(scores, 2, var), decreasing = TRUE)
  scores <- scores[, o]
  coef <- coef[o, ]
  
  # Remove if fewer than 5 genes
  scores.norm <- t(scores) / apply(scores, 2, mean)
  ranks_x <- t(apply(-t(scores.norm), 1, rank))
  ranks_y <- apply(-t(scores.norm), 2, rank)
  
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  
  modules <- apply(ranks_y, 2, function(m){
    a <- sort(m[is.finite(m)])
    a <- a[a == 1:length(a)]
    names(a)
  })
  
  l <- sapply(modules, length)
  keep <- (l >= 10)
  scores <- scores[, keep]
  coef <- coef[keep, ]
  
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