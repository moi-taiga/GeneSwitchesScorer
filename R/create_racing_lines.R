#' @title Identify the "racing lines"
#'
#' @description
#' Produces an estimate for the position on trajectory of each gene in each cell of a sample.
#'     This can be later aggregated to estimate the position of the sample along the trajectory.
#'
#' @param reduced_binary_counts_matrix a matrix of your samples binary gene expression.
#' @param gs_scorer_genes Switching genes which are evenly distributed through pseudotime.
#'
#' @return A list of matrices: A matrix for each cell where the columns represent progress through a trajectory,
#'     and the rows represent genes, values indicate a likely position of the cell upon the trajectory based that genes bianrized expression.
#' @export
#'
create_racing_lines <- function(reduced_binary_counts_matrix,gs_scorer_genes) {

  ## Input values:
  number_of_cells <- dim(reduced_binary_counts_matrix)[2]
  # Should we get the number of switching genes from here? or dim(reduced_binary_counts_matrix)[1]
  # It is possible for dim(reduced_binary_counts_matrix)[1] <  dim(gs_scorer_genes)[1]
  number_of_switching_genes <- dim(gs_scorer_genes)[1]


  # Different datasets store this information in different columns.. watch out.
  # Maybe this should be an input to the function?
  # as.numeric is needed/notneeded depending on dataset.
  switching_time <- as.numeric(gs_scorer_genes$switch_at_timeidx)
  switching_direction <- gs_scorer_genes$direction

  # Building the final list. (faster than building it dynamically.)
  all_patients_cells_scored <- vector("list", number_of_cells)
  names(all_patients_cells_scored) <- colnames(reduced_binary_counts_matrix)

  # Loop through all cells making matrices for each,
  #which represent likely position of a cell on a trajectory based on the expression of each gene.
  for (c in 1:number_of_cells) {
    # Build the matrix of 0's which has genes as rows and pseudotime indecies as columns.
   racing_mat <- matrix(0, nrow = number_of_switching_genes, ncol = 100)
   rownames(racing_mat) <- rownames(reduced_binary_counts_matrix)
    binarized_gene_expression_for_cell_c <- reduced_binary_counts_matrix[, c]

    up_indices <- which(binarized_gene_expression_for_cell_c == 1 & switching_direction == "up")
    down_indices <- which(binarized_gene_expression_for_cell_c == 1 & switching_direction == "down")
    not_up_indices <- which(binarized_gene_expression_for_cell_c == 0 & switching_direction == "up")
    not_down_indices <- which(binarized_gene_expression_for_cell_c == 0 & switching_direction == "down")

    if (length(up_indices) > 0) {
      racing_mat[up_indices, switching_time[up_indices]:100] <- 1
    }

    if (length(down_indices) > 0) {
      racing_mat[down_indices, 1:switching_time[down_indices]] <- 1
    }

    if (length(not_up_indices) > 0) {
      racing_mat[not_up_indices, 1:switching_time[not_up_indices]] <- 1
    }

    if (length(not_down_indices) > 0) {
      racing_mat[not_down_indices, switching_time[not_down_indices]:100] <- 1
    }

    all_patients_cells_scored[[c]] <- racing_mat
  }
  return(all_patients_cells_scored)
}
