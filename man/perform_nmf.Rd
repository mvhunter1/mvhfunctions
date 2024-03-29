% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_nmf.R
\name{perform_nmf}
\alias{perform_nmf}
\title{perform_nmf}
\usage{
perform_nmf(seurat_obj, rank = 20, cluster_name)
}
\arguments{
\item{seurat_obj}{Seurat object.}

\item{rank}{Rank for NMF, default = 20.}

\item{cluster_name}{Name of cluster or clusters you want to run NMF on.}
}
\value{
NMF results list with matrix, modules, scores and coef.
}
\description{
Perform NMF and create gene modules from expression matrix. Slightly modified code from Reuben and Dalia.
}
