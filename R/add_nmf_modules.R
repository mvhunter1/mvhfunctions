#' @title add_nmf_modules
#' @description Add NMF module scores to Seurat object metadata.
#' @param seurat_obj Seurat object.
#' @param nmf_results nmf results, ouput from perform_nmf.
#' @export
#' @return Seurat object.

add_nmf_modules <- function(seurat_obj, nmf_results) {
  require(Seurat)
  
  scores <- nmf_results$scores
  coef <- nmf_results$coef
  coef.z <- scale(coef)
  
  # Look at the expression of the top genes in each factor in the dataset
  factor.names <- paste0('NMF.', 1:dim(coef)[1])
  n_factors <- length(factor.names)
  
  for (ii in 1:dim(coef)[1]) {
    seurat_obj <- AddMetaData(seurat_obj, as.data.frame(coef[ii,]), col.name = factor.names[ii])
  }
  return(seurat_obj)
}