#' Title
#'
#' @param fibroblast_flat flattend sample data
#'
#' @return nice plot highlighting the probable position of your sample on your trajectory
#' @export
#'
gss_output_plot <- function(fibroblast_flat){
  #TODO: include an easy option to compare two samples.
  plot(x = 1:ncol(fibroblast_flat), y = fibroblast_flat[1,], type = "h", xlab = "Pseudotime Index", ylab = "Cell Position Likelyhood", main = "Trajectory Progress of \"fibroblast\"")
  abline(v = which.max(colSums(fib_flat)), lwd = 1, col = "Red")
}
