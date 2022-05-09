#' @title convert_seurat_obj_to_human
#' @description convert either Visium or scRNA-seq Seurat object from fish genes to human.
#' @param seurat_obj Seurat object.
#' @export
#' @return Seurat object with human genes instead of fish. Image data will be preserved if applying this function to Visium data.


convert_seurat_obj_to_human <- function(seurat_obj) {
  
  require(Seurat)
  require(tidyverse)
  require(data.table)
  
  message('Loading conversion table...')
  # load conversion table
  convert.table.Z11 <- tryCatch({
    read.delim("GRCz11_to_HS.txt")
  }, error = function(e) {
    stop('GRCz11_to_HS.txt must be in the working directory to use this function.')
  })
    
  # filter for DIOPT score higher than 6 (may have to go higher?)
  fish.human.convert.Z11 <- convert.table.Z11[convert.table.Z11$DIOPT_Score > 6, ] %>% select(Zebrafish_Symbol, Human_Symbol, Weighted_Score)
  
  # determine which slots are present in Seurat obj
  slots <- names(seurat_obj@assays)
  
  # determine if it is a Visium or single cell dataset.
  if ('Spatial' %in% slots) {
    dataset <- 'Visium'
    message('Visium dataset detected.')
  } else if ('RNA' %in% slots) {
    dataset <- "scRNA-seq"
    message('scRNA-seq dataset detected.')
  }
  
  if ('integrated' %in% slots) {
    message('Integrated slot detected...')
    message('Extracting integrated data...')
    # convert integrated data slot
    integrated_data <- GetAssayData(seurat_obj, assay = "integrated", slot = "data") %>%
      data.frame() %>%
      rownames_to_column(var = "fish_gene") %>%
      mutate(mean = rowMeans(.[,-1])) # add new column with mean expression
    
    message('Merging with conversion table...')
    # merge with conversion table
    integrated_data_merged <- merge(x = integrated_data %>% dplyr::select(fish_gene, mean),
                                    y = fish.human.convert.Z11,
                                    by.x = 'fish_gene',
                                    by.y = 'Zebrafish_Symbol')
    
    message('Finding and removing duplicated genes...')
    # find duplicated genes
    duplicated_human <- integrated_data_merged$Human_Symbol %>%
      .[duplicated(.)] %>%
      as.character() %>%
      unique() %>%
      as.list()
    
    # remove duplicates
    remove_dups_human <- lapply(duplicated_human, function(x) {
      integrated_data_merged[integrated_data_merged$Human_Symbol %in% x,] %>% # add associated expression data for each duplicated gene
        arrange(-mean) %>% # arrange from highest - lowest expression
        slice(1) # take the highest expressed gene
    })
    remove_dups_human <- rbindlist(remove_dups_human)
    
    # merge back with non duplicated genes
    integrated_data_human <- rbind(remove_dups_human,
                                   integrated_data_merged[!integrated_data_merged$Human_Symbol %in% unlist(duplicated_human),]) %>% # non duplicated genes 
      merge(x = .,
            y = integrated_data,
            by = c("mean", "fish_gene"))
    
    # Now extract SCT data for the same genes.
    message('Extracting SCT data...')
    SCT_data <- GetAssayData(seurat_obj, assay = "SCT", slot = "data") %>%
      data.frame() %>%
      rownames_to_column(var = "fish_gene")
    
    SCT_data_human <- merge(x = integrated_data_human %>% dplyr::select(fish_gene, Human_Symbol),
                            y = SCT_data,
                            by = 'fish_gene') %>%
      dplyr::select(-fish_gene) %>%
      column_to_rownames(var = "Human_Symbol")
    
    # Now extract counts data for the same genes.
    message('Extracting counts data...')
    if (dataset == 'Visium') {
      counts_data <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts") %>%
        data.frame() %>%
        rownames_to_column(var = "fish_gene")
    } else {
      counts_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts") %>%
        data.frame() %>%
        rownames_to_column(var = "fish_gene")
    }

    counts_data_human <- merge(x = integrated_data_human %>% dplyr::select(fish_gene, Human_Symbol),
                            y = counts_data,
                            by = 'fish_gene') %>%
      dplyr::select(-fish_gene) %>%
      column_to_rownames(var = "Human_Symbol")
    
    
    message('Creating new Seurat object...')
    # create new Seurat object with human genes
    integrated_data_human <- integrated_data_human %>%
      dplyr::select(-fish_gene, -mean, -Weighted_Score) %>%
      column_to_rownames(var = "Human_Symbol")
    
    # workaround to fix cell barcodes
    colnames(integrated_data_human) <- gsub(pattern = "\\.", replacement = "-", x = colnames(integrated_data_human))
    colnames(SCT_data_human) <- gsub(pattern = "\\.", replacement = "-", x = colnames(SCT_data_human))
    colnames(counts_data_human) <- gsub(pattern = "\\.", replacement = "-", x = colnames(counts_data_human))
    
   seurat_obj_human <- CreateSeuratObject(counts = integrated_data_human,
                                           assay = "integrated",
                                           meta.data = seurat_obj[[]]) # add metadata from previous Seurat obj
    
    #### IF THIS IS A VISIUM DATASET - NEED TO ADD THE EXPRESSION DATA IN A SLIGHTLY DIFFERENT WAY TO PRESERVE IMAGE DATA
    if (dataset == 'Visium') {
      message('Preserving image data...')
      # add data back to original Seurat object by overwriting the integrated slot. This is kinda messy but it works
      seurat_obj_human_2 <- seurat_obj
      integrated_data_to_add <- GetAssayData(seurat_obj_human[['integrated']], slot = 'data')
      seurat_obj_human_2[['integrated']] <- CreateAssayObject(data = integrated_data_to_add)
      seurat_obj_human_2[['SCT']] <- CreateAssayObject(data = as.matrix(SCT_data_human))
      seurat_obj_human_2[['Spatial']] <- CreateAssayObject(counts = as.matrix(counts_data_human))
      message('Done!')
      return(seurat_obj_human_2)
    } else {
      seurat_obj_human[['SCT']] <- CreateAssayObject(data = as.matrix(SCT_data_human))
      seurat_obj_human[['RNA']] <- CreateAssayObject(counts = as.matrix(counts_data_human))
      message('Adding PCA and UMAP embeddings...')
      # add PCA and UMAP loadings from previous Seurat obj
      seurat_obj_human[["umap"]] <- CreateDimReducObject(embeddings = seurat_obj@reductions$umap@cell.embeddings,
                                                         assay = "integrated")
      seurat_obj_human[["pca"]] <- CreateDimReducObject(embeddings = seurat_obj@reductions$pca@cell.embeddings,
                                                        assay = "integrated")
      message('Done!')
      return(seurat_obj)
    }
  }
}


  



  
 
  