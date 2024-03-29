% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_gene_names_deseq2.R
\name{add_gene_names_deseq2}
\alias{add_gene_names_deseq2}
\title{add_gene_names_deseq2}
\usage{
add_gene_names_deseq2(
  deseq_results = NULL,
  deseq_norm_matrix = NULL,
  counts_matrix = NULL
)
}
\arguments{
\item{deseq_results}{deseq2 results data frame.}

\item{deseq_norm_matrix}{normalized expression matrix.}

\item{counts_matrix}{counts matrix.}
}
\value{
named expression data.
}
\description{
Add gene names to deseq2 results, counts matrix or normalized expression matrix. Can only do one at a time for now.
}
