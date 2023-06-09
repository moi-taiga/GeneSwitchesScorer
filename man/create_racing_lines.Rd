% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_racing_lines.R
\name{create_racing_lines}
\alias{create_racing_lines}
\title{Identify the "racing lines"}
\usage{
create_racing_lines(reduced_binary_counts_matrix, gss_genes)
}
\arguments{
\item{reduced_binary_counts_matrix}{a matrix of your samples binary gene expression.}

\item{gss_genes}{Switching genes which are evenly distributed through pseudotime.}
}
\value{
A list of matrices: A matrix for each cell where the columns represent progress through a trajectory,
and the rows represent genes, values indicate a likely position of the cell upon the trajectory based that genes bianrized expression.
}
\description{
Produces an estimate for the position on trajectory of each gene in each cell of a sample.
This can be later aggregated to estimate the position of the sample along the trajectory.
}
