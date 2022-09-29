#' @title nice_bayesspace_feature_plot
#' @description Plot enhanced expression of a given gene resulting from BayesSpace analysis.
#' @param enhanced_obj Name of enhanced BayesSpace object.
#' @param normal_obj Optional: name of Seurat object to show comparison.
#' @param features Genes to plot expression of.
#' @param cols Colormap for plotting.
#' @param show_comparison if T, will plot both "normal" expression in Seurat object as well as enhanced expression.
#' @param show_title show titles on plots?
#' @param plot_enhanced if F, will only plot the Seurat expression.
#' @param n_col number of columns for final plot.
#' @export
#' @return plots.

nice_bayesspace_feature_plot <- function(enhanced_obj, normal_obj = NULL, features, cols = brewer.rdpu(n = 100), show_comparison = F, show_title = F, plot_enhanced = T, n_col = NULL) {
  
  require(BayesSpace)
  require(pals)
  require(cowplot)
  
  if (show_comparison) {
    if (is.null(normal_obj)) {
      stop('Visium object must be input to show comparison.')
    }
  }
  
  name_obj <- deparse(substitute(enhanced_obj))
  
  if (length(features) == 1) {
    if (show_comparison) {
      # plot normal resolution image
      p1 <- tryCatch({
        featurePlot(normal_obj, feature = features, size = 0) + scale_fill_gradientn(colours = cols)
      }, error = function(e) {
        stop(paste(features, "not found in BayesSpace object"))
      })
      # plot enhanced image
      p2 <- featurePlot(enhanced_obj, feature = features, size = 0) + scale_fill_gradientn(colours = cols)
      if (show_title) {
        p1 <- p1 + ggtitle("Visium") + 
          theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        p2 <- p2 + ggtitle("BayesSpace enhanced resolution") + 
          theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
      }
      plot <- plot_grid(p1,p2, ncol = 2)
    } else if (plot_enhanced == F) {
      plot <- tryCatch({
        featurePlot(normal_obj, feature = features, size = 0) + scale_fill_gradientn(colours = cols)
      }, error = function(e) {
        stop(paste(features, "not found in BayesSpace object"))
      })
    } else {
      plot <- tryCatch({
        featurePlot(enhanced_obj, feature = features, size = 0) + scale_fill_gradientn(colours = cols)
      }, error = function(e) {
        stop(paste(features, "not found in BayesSpace object"))
      })
    }
    return(plot)
    
  } else {
    enhanced_plots <- NULL
    
    if (show_comparison) {
      
      # first plot normal plots
      normal_plots <- NULL
      for (feature in features) {
        normal_plot <- tryCatch({
          featurePlot(normal_obj, feature = feature, size = 0) + scale_fill_gradientn(colours = cols)
        }, error = function(e) {
          message(paste(feature, "not found in BayesSpace object"))
        })
        if (inherits(normal_plot, "error")) {
          next
        }
        else {
          normal_plots[[feature]] <- normal_plot
        }
      }
      
      # then plot enhanced plots
      for (feature in features) {
        enhanced_plot <- tryCatch({
          featurePlot(enhanced_obj, feature = feature, size = 0) + scale_fill_gradientn(colours = cols)
        }, error = function(e) {
          # message(paste(feature, "not found in BayesSpace object"))
        })
        if (inherits(enhanced_plot, "error")) {
          next
        }
        else {
          enhanced_plots[[feature]] <- enhanced_plot
        }
      }
      
      normal_plots <- plot_grid(plotlist = normal_plots, ncol = length(normal_plots))
      enhanced_plots <- plot_grid(plotlist = enhanced_plots, ncol = length(enhanced_plots))
      plots <- plot_grid(normal_plots, enhanced_plots, nrow = 2)
      return(plots)
      
    } else {
      
      for (feature in features) {
        enhanced_plot <- tryCatch({
          featurePlot(enhanced_obj, feature = feature, size = 0) + scale_fill_gradientn(colours = cols)
        }, error = function(e) {
          message(paste(feature, "not found in BayesSpace object"))
        })
        if (inherits(enhanced_plot, "error")) {
          next
        }
        else {
          enhanced_plots[[feature]] <- enhanced_plot
        }
      }
      
      if (is.null(n_col)) {
        plots <- plot_grid(plotlist = enhanced_plots)
      } else {
        plots <- plot_grid(plotlist = enhanced_plots, ncol = n_col)
      }
      
      return(plots)
    }
  }
}
