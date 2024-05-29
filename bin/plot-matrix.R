#!/usr/bin/env Rscript

###############################################################################
# -*- encoding: UTF-8 -*-                                                     #
# Author: Jesse C. Chen (jessekelighine.com)                                  #
# Description: Plot Matrix                                                    #
#                                                                             #
###############################################################################
require(ggplot2)
###############################################################################

plot.matrix <- function ( mat, mutual.highlight=FALSE, over.half.highlight=FALSE ) {
  output <- data.frame(
    row = ( 1:length(mat) - 1 ) %%  nrow(mat) + 1,
    col = ( 1:length(mat) - 1 ) %/% nrow(mat) + 1,
    value = as.vector(as.numeric(mat)),
    mutual = as.vector(mat==t(mat) & ( mat!=0 | t(mat)!=0 ))
  ) |>
    transform(over.half = value >= 0.5) |> 
    ggplot() +
    geom_raster(aes(x=col, y=row, fill=value)) +
    scale_y_reverse(breaks = 1:nrow(mat)) +
    scale_x_continuous(breaks = 1:nrow(mat)) +
    scale_fill_gradient(low = "white", high = "red")
  if ( mutual.highlight ) {
    output <- output +
    geom_point(aes(x=col, y=row, color=mutual), size=1) +
    scale_color_manual(values = c("transparent", "white"))
  }
  if ( over.half.highlight ) {
    output <- output +
    geom_point(aes(x=col, y=row, color=over.half), size=1) +
    scale_color_manual(values = c("transparent", "black"))
  }
  output +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
}

### Testing ###################################################################

if ( sys.nframe() == 0 ) {
  size <- 10
  sample(0:1, size=size^2, replace=TRUE) |> 
    (`+`)( rnorm(size^2) ) |> 
    matrix(nrow=size) |> 
    plot.matrix(mutual.highlight=FALSE)
}
