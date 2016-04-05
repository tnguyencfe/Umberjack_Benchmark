library(grid) # for ggplot in grid formation  (unit.pmax)
library(gridExtra)  # for arranging ggplot in grid formation
library(gtable)  # for arranging in grid formation
library(ggplot2)

# Stolen from http://stackoverflow.com/questions/15016995/how-to-align-multiple-ggplot2-plots-and-add-shadows-over-all-of-them
# http://stackoverflow.com/questions/26159495/align-multiple-ggplot-graphs-with-and-without-legends
#+ fig.width=20, fig.height=12

AlignPlots <- function(..., samePlotSize=FALSE) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  max.heights <- do.call(unit.pmax, lapply(plots.grobs, "[[", "heights"))

  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
#   if (samePlotSize) {  # same heights for all plots
#     plots.grobs.eq.widths <- lapply(plots.grobs, function(g){
#       panels <- g$layout$t[grep("panel", g$layout$name)]
#       g$heights[panels] <- max.heights[panels]
#       
#       background <- as.vector(t(g$layout[grep("background", g$layout$name), c("t", "l", "b", "r")]))
#       g$heights[background] <- max.heights[background]
#       
#       return (g)
#     })
#   }
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}



# This removes outliers so that we can visualize the majority of points 
outlier_range <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE, names=FALSE)
  fudge <- 1.5 * IQR(x, na.rm = TRUE)
  return (c(lower=max(qnt[1] - fudge, min(x, na.rm=TRUE)), 
            upper=min(qnt[2] + fudge, max(x, na.rm=TRUE))))
}
