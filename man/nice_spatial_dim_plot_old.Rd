% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nice_spatial_dim_plot_old.R
\name{nice_spatial_dim_plot_old}
\alias{nice_spatial_dim_plot_old}
\title{nice_spatial_dim_plot_old}
\usage{
nice_spatial_dim_plot_old(
  seurat_obj,
  group.by,
  im_alpha = 0,
  pt.size = 1.4,
  stroke = 0,
  cols = NULL,
  label = F,
  show_legend = T
)
}
\arguments{
\item{seurat_obj}{Seurat object.}

\item{group.by}{what to colour the points by, usually a column in the Seurat object metadata.}

\item{im_alpha}{set to 1 to plot the tissue image, 0 otherwise.}

\item{pt.size}{size of plotted points on spatial array.}

\item{stroke}{linewidth to outline plotted points in black (default = no outline).}

\item{cols}{optional: vector of colours for each plotted group.}

\item{label}{if T, will label groups with text on plot.}

\item{show_legend}{if T, will include the legend but won't label points w/text on plot.}
}
\value{
SpatialPlot.
}
\description{
nicer looking version of the Seurat function SpatialDimPlot.
}
