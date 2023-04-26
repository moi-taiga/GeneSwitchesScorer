#' @title Identify the "racing lines"
#'
#' @description
#' Produces an estimate for the position on trajectory of each gene in each cell of a sample.
#'     This can be later aggregated to estimate the position of the sample along the trajectory.
#'
#' @param reduced_binary_counts_matrix a matrix of your samples binary gene expression.
#' @param gs_scorer_genes Switching genes which are evenly distributed through pseudotime
#'
#' @return
#' @export
#'
#' @examples create_racing_lines(reduced_binary_counts_matrix,gs_scorer_genes)
#'
create_racing_lines <- function(reduced_binary_counts_matrix,gs_scorer_genes) {
  all_patients_cells_scored<-list()
  number_of_cells<-dim(reduced_binary_counts_matrix)[2]
  number_of_switching_genes<-dim(gs_scorer_genes) [1]
  #for each patient cell C
  for (c in 1:number_of_cells){
    #print(paste(c,"out of", number_of_cells, "cells"))
    racing_mat<-matrix(0,nrow = number_of_switching_genes, ncol = 100)

    #for each switching gene G
    for (g in 1:number_of_switching_genes){
      #print(paste(g,"out of", number_of_switching_genes, "switching genes"))
      #find out the switch time Gt, and direction Gd
      switching_time <- as.numeric(gs_scorer_genes[g,3])
      switching_direction <- gs_scorer_genes[g,2]
      #find out if its expressed Ct
      is_expressed <- reduced_binary_counts_matrix[g,c]
      #If Ct = TRUE
      if(is_expressed == 1){
        #If Gd = UP
        if(switching_direction == "up"){
          # [Gt:100] = 1
          racing_mat[g,switching_time:100]<-1
        }else{
          #If Gd = DOWN
          # [0:Gt] = 1
          racing_mat[g,0:switching_time]<-1
        }
        #If Ct = FALSE
      }else{
        if(switching_direction == "up"){
          #If Gd = UP
          #[0:Gt] = 1
          racing_mat[g,0:switching_time]<-1
        }else{
          #If Gd = DOWN
          # [Gt:100] = 1
          racing_mat[g,switching_time:100]<-1
        }
      }
    }
    all_patients_cells_scored[[c]]<-racing_mat
  }
  return(all_patients_cells_scored)
}
