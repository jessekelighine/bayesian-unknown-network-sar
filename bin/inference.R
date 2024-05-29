#!/usr/bin/env Rscript

###############################################################################
# -*- encoding: UTF-8 -*-                                                     #
# Author: Jesse C. Chen (jessekelighine.com)                                  #
# Description: Posterior Inference                                            #
#                                                                             #
# Last Modified: 2024-05-06                                                   #
###############################################################################
options(scipen=999999)
library(txtplot)
###############################################################################

rm(list=ls(all.name=TRUE))

input <- commandArgs(trailingOnly=TRUE) # Read from STDIN
TRIAL <- input[1]
plot.matrix.size <- if ( length(input) < 2 ) 12 else as.numeric(input[2])

###############################################################################

read.dir <- paste0("output-param-", TRIAL)
chain.length <- list.files(read.dir) |>
  Filter( f = \ (.) grep(., pattern="\\.rds$") ) |>
  length()
usable.length <- floor(min(3e3, chain.length) * 0.9)
c(chain.length=chain.length, usable.length=usable.length)

### Read Data #################################################################

name.rds <- function ( iter ) paste0("param-", iter, ".rds")
output.chain <- tail(1:chain.length, n=usable.length) |>
  lapply( \ (iter) file.path(read.dir, name.rds(iter)) ) |>
  lapply( readRDS )

### Txtplots ##################################################################

output.chain |>
  # lapply( `[[`, "beta" ) |> lapply( `[`, 4 ) |>
  lapply( `[[`, "lambda" ) |>
  unlist() |>
  txtplot(x=1:usable.length, y=_)
  # txtdensity()

### Adjacency Matrix ##########################################################

source("plot-matrix.R")

posterior.W <- output.chain |>
  lapply( `[[`, "net" ) |>
  Reduce(f=`+`) |>
  (`/`)( length(output.chain) )
posterior.W.plot <- plot.matrix(posterior.W, over.half.highlight=TRUE)

true.W <- file.path("data", paste0(TRIAL, ".rds")) |> readRDS() |> (`$`)("net")
# true.W.plot <- true.W |> plot.matrix(mutual.highlight=TRUE)

cat("MSE:", (true.W - posterior.W)^2 |> sum(), fill=TRUE)

ggsave(
  filename = file.path("output-figures", paste0(TRIAL, "-posterior-W.pdf")),
  plot     = posterior.W.plot,
  width    = plot.matrix.size,
  height   = plot.matrix.size,
  unit     = "cm"
)

### Regression Table ##########################################################

regression.table <- list()

regression.table$to.data.frame <- function ( output.chain ) {
  params.to.keep <- c("lambda","beta","sigma2","theta")
  beta.names  <- paste0("beta",  1:length(output.chain[[1]]$beta))
  theta.names <- paste0("theta", 1:length(output.chain[[1]]$theta))
  output.chain |>
    lapply( \ (param) param[names(param) %in% params.to.keep] ) |>
    lapply( \ (param) c(param, as.list(as.vector(param$beta))  |> (`names<-`)(beta.names)) ) |>
    lapply( \ (param) c(param, as.list(as.vector(param$theta)) |> (`names<-`)(theta.names)) ) |>
    lapply( `[[<-`, "beta", NULL ) |>
    lapply( `[[<-`, "theta", NULL ) |>
    lapply( as.data.frame ) |>
    do.call( rbind, args=_ )
}

regression.table$generate <- function (
  output.chain,
  qt.probs = c(0.025, 0.975)
) {
  output.chain |>
    regression.table$to.data.frame() |>
    as.list() |>
    lapply( \ (param) {
      list(
        post.mean = mean(param),
        post.qt = quantile(param, probs=qt.probs)
      )
    } ) |>
    lapply( unlist ) |>
    lapply( sapply, sprintf, fmt="%1.4f" ) |>
    lapply( \ (vec) {
      new.names <- names(vec) |> sapply(gsub, pattern="%", replacement="")
      names(vec) <- new.names
      vec
    } ) |>
    lapply( t ) |>
    lapply( as.data.frame ) |>
    do.call( rbind, args=_ )
}

regression.table$generate(output.chain) |>
  capture.output(file=file.path("output-tables", paste0(TRIAL,".txt")))

quit()

### Tidyverse: Posterior Distribution #########################################

library(conflicted)
library(tidyverse)

output.chain |>
  lapply( \ (param) c(param$lambda, param$sigma2) ) |>
  unlist() |>
  matrix(ncol=2, byrow=TRUE) |>
  as_tibble() |>
  ggplot() +
  aes(x=V1, y=V2) +
  geom_point()

pdf("Rplots.pdf")
output.chain |>
  head(501) |> tail(1) |> 
  (`[[`)(1) |> 
  (`$`)(net) |> 
  plot.matrix()
dev.off()
