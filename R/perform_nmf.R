#' @title perform_nmf
#' @description Perform NMF and create gene modules from expression matrix. Slightly modified code from Reuben and Dalia.
#' @param expression_matrix Expression matrix, commonly variable genes, for cell type of interest.
#' @param rank Rank for NMF, default = 20. Can be a vector of possible ranks but can only use multiple ranks if use_fast_method = F.
#' @param save_NMF_results Save NMF results matrix in working dir?
#' @param filename Filename to save NMF results as.
#' @param use_fast_method If T, will use the fast NMF method implemented in RcppML. Can only be used for a single rank.
#' @export
#' @return NMF results list with matrix, modules, scores and coef.

perform_nmf <- function(expression_matrix, rank = 20, save_NMF_results = F, filename = NULL, use_fast_method = F) {
  
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
  if (use_fast_method) {
    
    if (length(ranks) > 1) {
      stop('Multiple ranks cannot be used with use_fast_method.')
    }
    
    require(RcppML) 
    require(Matrix)
    nmf_results <- RcppML::nmf(A = Matrix::Matrix(expression_matrix_filt, sparse = T), 
                               k = rank, 
                               seed = 'ica') 
  } else {
    require(NMF)
    require(fastICA)
    message('Using slow NMF method. Warning: this may take hours to run...')
    nmf_results <- NMF::nmf(x = expression_matrix_filt, 
                            rank = rank, 
                            seed = 'ica', 
                            #.options = list(parallel=T),
                            method = 'nsNMF')
  }
  t2 <- Sys.time() 
  message('NMF calculations complete!')
  difftime(t2,t1)
  
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
  if (use_fast_method) {
    scores <- nmf_results$w
    coef <- nmf_results$h
    rownames(scores) <- rownames(expression_matrix_filt) # critical!!!!
    colnames(coef)<- colnames(expression_matrix_filt) # critical!!!!
  } else {
    scores <- basis(nmf_results)
    coef <- coefficients(nmf_results)
  }
  
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