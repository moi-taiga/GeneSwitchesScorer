% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_gene_expression_for_switching_genes.R
\name{filter_gene_expression_for_switching_genes}
\alias{filter_gene_expression_for_switching_genes}
\title{Filter Gene Expression for Switching Genes}
\usage{
filter_gene_expression_for_switching_genes(
  binary_counts_matrix,
  gs_scorer_genes
)
}
\arguments{
\item{binary_counts_matrix}{a binary expression matrix from your sample.}

\item{gs_scorer_genes}{Switching genes which are evenly distributed through pseudotime..}
}
\value{
a reduced binary expression matrix filtered to only include selected switching genes
}
\description{
Create a reduced binary expression matrix for only the selected switching genes,
binary_counts_matrix is from the sample DATA and gs_scorer_genes is from Atlas Data.
}
