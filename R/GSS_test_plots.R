library(ggrepel)


# To look at the racing lines of only one cell you need to make sure it still has the right matrix format.
fib_red <- matrix(fibroblast_reduced[, 1], ncol = 1, dimnames = list(rownames(fibroblast_reduced), colnames(fibroblast_reduced)[1]))
fib_lines   <- create_racing_lines(fib_red, gss_genes)
#
plot_timeline_ggplot(gss_genes, timedata = colData(reference_glm.gs)$Pseudotime, txtsize = 3)

###################################### Use Pdeudotime IDX ############################################################################################

tml <- gss_genes
fib_red <- fibroblast_reduced
#racing_lines <- function(tml, iffulltml = TRUE, txtsize = 3.5) {

  # Convert tml to a data frame
  tml <- as.data.frame(tml)

  # Order tml by switch_at_timeidx column
  tml <- tml[order(tml$switch_at_timeidx), ]

  # Add a new column direction_num and set it to -1
  tml$direction_num <- -1

  # If "up" is present in the direction column, set direction_num to 1
  if ("up" %in% tml$direction) {
    tml[tml$direction == "up", ]$direction_num <- 1
  }

  # Generate pseudotime data based on the specified parameters
  if (iffulltml) {
    timeidx_df <- data.frame(timeidx_range = c(0, 25, 50, 75, 100), timeidx_labs = c(0, 25, 50, 75, 100))
  } else {
    timeidx_step <- (max(tml$switch_at_timeidx) - min(tml$switch_at_timeidx))/4
    timeidx_range <- seq(min(tml$switch_at_timeidx), max(tml$switch_at_timeidx), by = timeidx_step)
    timeidx_df <- data.frame(timeidx_range, timeidx_labs = round(timeidx_range, 1))
  }

  # Create the initial ggplot object with x and y aesthetics, color, and labels
  tml_plot <- ggplot(tml, aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(tml))) +
    geom_point(size = txtsize/3) + xlab("Time-Index") + ylab("Quality of fitting (R^2)")

  # Add the classic theme to the plot
  tml_plot <- tml_plot + theme_classic()

  # Add a horizontal black line for the timeline
  tml_plot <- tml_plot + geom_hline(yintercept = 0, color = "black", linewidth = 0.6)

  # Add labels for pseudotime on the plot
  tml_plot <- tml_plot + geom_label(data = timeidx_df, aes(x = timeidx_range, y = 0, label = timeidx_labs),
                                    size = (txtsize-0.5), color = "black")

  # Add text labels with repulsion to avoid overlap
  tml_plot <- tml_plot + geom_text_repel(aes(x = switch_at_timeidx, y = pseudoR2s * direction_num, label = rownames(tml)),
                                         size = txtsize, show.legend = FALSE)

  # Customize the theme and legend appearance
  tml_plot <- tml_plot + theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(10, "pt"),
                               text = element_text(size = 12, family = "Helvetica"))

  # title the plot the name of the cell.
  tml_plot <- tml_plot + ggtitle()

  #
  tml_plot
  # if G1 is expressed in C1 and the switch is UP then draw the line to the right.
  # if G1 is expressed in C1 and the switch is Down then draw the line to the Left.
  # if G1 is NOT expressed in C1 and the switch is UP then draw the line to the right.
  # if G1 is NOT expressed in C1 and the switch is Down then draw the line to the Left.

  #this needs to be changed to refer to the sample.

  lined_plot <- tml_plot

  for (i in 1:dim(fib_red)[1]) {
    if ((fib_red[rownames(tml)[i], 1] == 0) && (tml$direction[i] == "up")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = 0, xend = tml$switch_at_timeidx[i], y = tml$pseudoR2s[i], yend = tml$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }

    if ((fib_red[rownames(tml)[i], 1] == 1) && (tml$direction[i] == "up")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = tml$switch_at_timeidx[i], xend = 100, y = tml$pseudoR2s[i], yend = tml$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }

    if ((fib_red[rownames(tml)[i], 1] == 0) && (tml$direction[i] == "down")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = tml$switch_at_timeidx[i], xend = 100, y = -tml$pseudoR2s[i], yend = -tml$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }

    if ((fib_red[rownames(tml)[i], 1] == 1) && (tml$direction[i] == "down")) {
      lined_plot <- lined_plot + geom_segment(aes_string(x = 0, xend = tml$switch_at_timeidx[i], y = -tml$pseudoR2s[i], yend = -tml$pseudoR2s[i]),
                                              color = "blue", size = 0.8)
    }
  }

  lined_plot




  # Return the final plot
  return(lined_plot)
}
