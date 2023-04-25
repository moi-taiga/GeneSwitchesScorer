#' Title flatten_cell_frequencies_moi
#'
#' @description
#' Combining all cells racing lines with out binerization (Moi's way).
#'
#'
#' @param list_of_cell_position_frequencies
#'
#' @return
#' @export
#'
#' @examples
flatten_cell_frequencies_moi <- function(list_of_cell_position_frequencies) {
  #making an empty flat matrix
  all_patient_cells_scored_flat <- matrix(0, nrow = 1, ncol = 100)
  #
  length_of_list <- length(list_of_cell_position_frequencies)
  # for every cells matrix calculate the colsums and add them to the flat matrix
  for (i in 1:length_of_list) {
    all_patient_cells_scored_flat <- all_patient_cells_scored_flat + colSums(list_of_cell_position_frequencies[[i]])
  }
  # divide the values in the flat matrix by the number of cells in an attempt to make the yaxis more informative.
  # This may be the wrong approach but it shouldnt change the shape of the final plot
  all_patient_cells_scored_flat <- all_patient_cells_scored_flat/length_of_list

  return(all_patient_cells_scored_flat)
}
