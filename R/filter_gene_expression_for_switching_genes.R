#' Filter Gene Expression for Switching Genes
#'
#' @description
#' Create a reduced binary expression matrix for only the selected switching genes,
#'     binary_counts_matrix is from the Patient DATA and gs_scorer_genes is from Atlas Data.
#' @param binary_counts_matrix
#' @param gs_scorer_genes
#'
#' @return
#' @export
#'
#' @examples
#'
filter_gene_expression_for_switching_genes<-function(binary_counts_matrix, gs_scorer_genes) {
  indices_of_switching_genes<-which(rownames(binary_counts_matrix) %in% gs_scorer_genes[,1])
  reduced_binary_counts_matrix<- binary_counts_matrix[indices_of_switching_genes,]
  gs_scorer_genes_to_keep<-which(gs_scorer_genes[,1] %in% rownames(reduced_binary_counts_matrix))
  gs_scorer_genes<- gs_scorer_genes[gs_scorer_genes_to_keep,]
  #print(gs_scorer_genes_to_keep)
  returnlist<-list(reduced_binary_counts_matrix, gs_scorer_genes)
  return(returnlist)
}
