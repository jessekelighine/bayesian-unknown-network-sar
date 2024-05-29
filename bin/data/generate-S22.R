#!/usr/bin/env Rscript

###############################################################################
# -*- encoding: UTF-8 -*-                                                     #
# Author: Jesse C. Chen (jessekelighine.com)                                  #
# Description: Generate Mock Data                                             #
#                                                                             #
###############################################################################
options(scipen=999999)
set.seed(2024)
###############################################################################
# library(ergm)
###############################################################################

rm(list=ls(all.name=TRUE)); gc()

na.replace <- function ( object, with=0 ) ifelse(is.na(object), with, object)
row.normalize <- function ( matrix ) na.replace(matrix / rowSums(matrix))

###############################################################################

TRIAL <- "S22"

data <- list()
data$lambda <- 0.1
data$beta1 <- c(2, 0.5)
data$beta2 <- c(2, 0.5)
data$theta <- c(edges=-1.8, mutual=2.3, ctriad=-0.1)
data$beta <- c(data$beta1, data$beta2)
data$sigma2 <- 2
data$panel_length <- 8
data$size <- 18

data$net <- readRDS("sampson.rds") |> as.matrix()
data$W <- row.normalize(data$net)

data$U <- 1:data$panel_length |> lapply( \ (time) rnorm(n=data$size, mean=0, sd=sqrt(data$sigma2)) )

.X2 <- sample(0:1, size=data$size, replace=TRUE)
data$X1 <- 1:data$panel_length |>
  lapply( \ (time) c(
    rnorm(n=data$size, mean=10, sd=10)
  ) ) |>
  lapply( matrix, ncol = 1 ) |>
  lapply( \ (matrix) cbind(matrix, .X2) )

data$X2 <- data$X1
# .X2 <- sample(0:1, size=data$size, replace=TRUE)
# data$X2 <- 1:data$panel_length |>
#   lapply( \ (time) c(
#     rnorm(n=data$size, mean=10, sd=10)
#   ) ) |>
#   lapply( matrix, ncol = 1 ) |>
#   lapply( \ (matrix) cbind(matrix, .X2) )

data$Y <- 1:data$panel_length |>
  lapply( function (time) {
    M <- solve(diag(data$size) - data$lambda * data$W)
    H <- cbind(data$X1[[time]], data$W %*% data$X2[[time]])
    M %*% ( H %*% data$beta + data$U[[time]] )
  } )

saveRDS(data, paste0(TRIAL, ".rds"))

### Testing ###################################################################

if ( sys.nframe() == 0 ) {
  source("../plot-matrix.R")
  mutual.net <- data$net==t(data$net) & (data$net==1 | t(data$net==1))
  plot.matrix(data$net + mutual.net)
  plot(as.network(data$net))
}
