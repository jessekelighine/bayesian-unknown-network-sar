#!/usr/bin/env Rscript

###############################################################################
# -*- encoding: UTF-8 -*-                                                     #
# Author: Jesse C. Chen (jessekelighine.com)                                  #
# Description:                                                                #
#                                                                             #
# Last Modified: OOOO-OO-OO                                                   #
###############################################################################
options(scipen=999999)
# options(datatable.verbose=FALSE,datatable.quiet=TRUE)
# options(tidyverse.quiet=TRUE)
set.seed(2024)
###############################################################################
# library(conflicted)
# library(data.table)
# library(tidyverse)
# library(network)
# library(Bergm)
###############################################################################

rm(list=ls(all.name=TRUE)); gc()

TRIAL <- "K11"

na.replace <- function ( object, with=0 ) ifelse(is.na(object), with, object)
row.normalize <- function ( matrix ) na.replace(matrix / rowSums(matrix))

###############################################################################

data <- list()
data$lambda <- 0.05
data$beta1 <- c(2, 1)
data$beta2 <- 0.5
data$theta <- c(
  edges = -3,
  kstar2 = 0.2,
  kstar3 = -0.1,
  triangle = 0.3
)
data$beta <- c(data$beta1, data$beta2)
data$sigma2 <- 5
data$panel_length <- 12
data$net <- readRDS("zach.rds") |> as.matrix()
data$size <- nrow(data$net)
data$W <- row.normalize(data$net)
data$U <- 1:data$panel_length |> lapply( \ (time) rnorm(n=data$size, mean=0, sd=sqrt(data$sigma2)) )
data$X1 <- 1:data$panel_length |>
  lapply( \ (time) {
    c(
      rnorm(n=data$size, mean=80, sd=15),
      sample(0:1, size=data$size, replace=TRUE)
    )
  } ) |>
  lapply( matrix, ncol=2 )
data$X2 <- data$X1 |> lapply(`[`,,1) |> lapply(matrix)
data$Y <- 1:data$panel_length |>
  lapply( function (time) {
    M <- diag(data$size) - data$lambda * data$W
    H <- cbind(data$X1[[time]], data$W %*% data$X2[[time]])
    solve(M) %*% ( H %*% data$beta + data$U[[time]] )
  } )

saveRDS(data, paste0(TRIAL, ".rds"))

### Testing ###################################################################

if ( sys.nframe() == 0 ) {

  panel.data <- readRDS("K1.rds")
  panel.data$beta

  source("../plot-matrix.R")
  mutual.net <- data$net==t(data$net) & (data$net==1 | t(data$net==1))
  plot.matrix(data$net + mutual.net)
  plot(as.network(data$net))

  pdf(file="zach.pdf")
  readRDS("zach.rds") |> as.network(directed=FALSE) |> plot()
  dev.off()


}
