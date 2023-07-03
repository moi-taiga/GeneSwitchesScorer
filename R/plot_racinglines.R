#' @title Identify and Visualise each cells position.
#'
#' @description
#' Produces a plot for each cell which helps visualise how GSS is predicting the cells postion.
#'
#' @param reduced_binary_counts_matrix a matrix of your samples binary gene expression.
#' @param gs_scorer_genes Switching genes which are evenly distributed through pseudotime.
#'
#' @return A list of matrices: A matrix for each cell where the columns represent progress through a trajectory,
#'     and the rows represent genes, values indicate a likely position of the cell upon the trajectory based that genes bianrized expression.
#' @import ggplot2
#' @import ggrepel
#'
#'
#' @export
#'



# To look at the racing lines of only one cell you need to make sure it still has the right matrix format.
fib_red <- matrix(fibroblast_reduced[, 1], ncol = 1, dimnames = list(rownames(fibroblast_reduced), colnames(fibroblast_reduced)[1]))
fib_lines   <- create_racing_lines(fib_red, gss_genes)
#


###################################### Use Pdeudotime IDX ############################################################################################

gss_genes <- gss_genes_tmp
reduced_binary_counts_matrix <- fibroblast_reduced
racinglines_timeline <- function(gss_genes, reduced_binary_counts_matrix, full_time_IDX = TRUE) {

  # Convert tml to a data frame
  gss_genes <- as.data.frame(gss_genes)

  # Order gss_genes by switch_at_timeidx column
  gss_genes <- gss_genes[order(gss_genes$switch_at_timeidx), ]

  # Add a new column direction_num and set it to -1
  gss_genes$direction_num <- -1

  # If "up" is present in the direction column, set direction_num to 1
  if ("up" %in% gss_genes$direction) {
    gss_genes[gss_genes$direction == "up", ]$direction_num <- 1
  }

  # Generate pseudotime data based on the specified parameters
  if (full_time_IDX) {
    timeidx_df <- data.frame(timeidx_range = c(0, 25, 50, 75, 100), timeidx_labs = c(0, 25, 50, 75, 100))
  } else {
    timeidx_step <- (max(gss_genes$switch_at_timeidx) - min(gss_genes$switch_at_timeidx))/4
    timeidx_range <- seq(min(gss_genes$switch_at_timeidx), max(gss_genes$switch_at_timeidx), by = timeidx_step)
    timeidx_df <- data.frame(timeidx_range, timeidx_labs = round(timeidx_range, 1))
  }

  # Create the initial ggplot object with x and y aesthetics, color, and labels
  timeline_plot <- ggplot(gss_genes, aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(gss_genes))) +
    geom_point(size = 1) + xlab("Time-Index") + ylab("Quality of fitting (R^2)")

  # Add the classic theme to the plot
  timeline_plot <- timeline_plot + theme_classic()

  # Add a horizontal black line for the timeline
  timeline_plot <- timeline_plot + geom_hline(yintercept = 0, color = "black", linewidth = 0.6)

  # Add labels for pseudotime on the plot
  timeline_plot <- timeline_plot + geom_label(data = timeidx_df, aes(x = timeidx_range, y = 0, label = timeidx_labs),
                                    size = (3), color = "black")

  # Add text labels with repulsion to avoid overlap
  timeline_plot <- timeline_plot + geom_text_repel(aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(gss_genes)),
                                         size = 3, show.legend = FALSE)

  # Customize the theme and legend appearance
  timeline_plot <- timeline_plot + theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(10, "pt"),
                               text = element_text(size = 12, family = "Helvetica"))

  # title the plot the name of the cell.
  timeline_plot <- timeline_plot + ggtitle(colnames(fib_red)[1])

  #
  timeline_plot
  # if G1 is expressed in C1 and the switch is UP then draw the line to the right.
  # if G1 is expressed in C1 and the switch is Down then draw the line to the Left.
  # if G1 is NOT expressed in C1 and the switch is UP then draw the line to the right.
  # if G1 is NOT expressed in C1 and the switch is Down then draw the line to the Left.

  #this needs to be changed to refer to the sample.

  lined_plot <- timeline_plot

  for (i in 1:dim(fib_red)[1]) {
    if ((fib_red[rownames(gss_genes)[i], 1] == 0) && (gss_genes$direction[i] == "up")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = 0, xend = gss_genes$switch_at_timeidx[i], y = gss_genes$pseudoR2s[i], yend = gss_genes$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }

    if ((fib_red[rownames(gss_genes)[i], 1] == 1) && (gss_genes$direction[i] == "up")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = gss_genes$switch_at_timeidx[i], xend = 100, y = gss_genes$pseudoR2s[i], yend = gss_genes$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }

    if ((fib_red[rownames(gss_genes)[i], 1] == 0) && (gss_genes$direction[i] == "down")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = gss_genes$switch_at_timeidx[i], xend = 100, y = -gss_genes$pseudoR2s[i], yend = -gss_genes$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }

    if ((fib_red[rownames(gss_genes)[i], 1] == 1) && (gss_genes$direction[i] == "down")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = 0, xend = gss_genes$switch_at_timeidx[i], y = -gss_genes$pseudoR2s[i], yend = -gss_genes$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }
  }

  lined_plot




  # Return the final plot
  return(lined_plot)
}
