#!/usr/bin/env Rscript

###############################################################################
# -*- encoding: UTF-8 -*-                                                     #
# Author: Jesse C. Chen (jessekelighine.com)                                  #
# Description: Simulation                                                     #
#                                                                             #
# Last Modified: 2024-01-18                                                   #
###############################################################################
options(scipen=999999)
set.seed(2026)
###############################################################################
library(Rcpp)
###############################################################################

### Timer #####################################################################

timer <- list()
timer$now <- Sys.time()
timer$tic <- function ( ... ) { cat(..., fill=...length()!=0); timer$now <<- Sys.time() }
timer$toc <- function ( ... ) { cat(..., format(difftime(Sys.time(),timer$now)), fill=TRUE) }

### Read Data #################################################################

TRIAL <- "S22"
panel.data <- file.path("data", paste0(TRIAL, ".rds")) |> readRDS()

data <- list()
data$size <- panel.data$size
data$panel_length <- panel.data$panel_length
data$X1 <- panel.data$X1
data$X2 <- panel.data$X2
data$Y <- panel.data$Y

### Helper Functions ##########################################################

sourceCpp(file.path("src", "simulation-helper.cpp"))
sourceCpp(file.path("src", paste0("simulation-", TRIAL, "-loglikelihood.cpp")))

### Initilizations ############################################################

init.param <- function ( data ) {
  output <- list()
  output$lambda <- 0
  output$beta <- rep(0, ncol(data$X1[[1]]) + ncol(data$X2[[1]]))
  output$net <- sample(0:1, size=data$size^2, replace=TRUE, prob=c(0.9,0.1)) |>
    matrix(nrow=data$size) |>
    (`diag<-`)(0)
  output$W <- RowNormalize(output$net)
  output$theta <- rep(0, length(panel.data$theta))
  output$sigma2 <- 1
  output
}

### Conditional Posterior Samplers ############################################

sample.posterior <- list()

sample.posterior$lambda <- function ( data, param ) {
  log.density <- function ( lambda ) LogLikelihoodAll(
    data = data,
    lambda = lambda,
    W = param$W,
    beta1 = param$beta[1:2],
    beta2 = param$beta[3:4],
    sigma2 = param$sigma2
  )
  propose <- function ( state ) runif(1) * 2 - 1
  propose.log.density <- function ( here, given ) 0
  chain.length <- sample(50:100, size=1)
  lambda.state <- rbeta(1, shape1=4, shape2=4) * 2 - 1
  for ( iteration in 1:chain.length ) {
    lambda.proposal <- propose(lambda.state)
    metropolis.ratio <- 0
    metropolis.ratio <- metropolis.ratio + log.density(lambda.proposal)
    metropolis.ratio <- metropolis.ratio - log.density(lambda.state)
    metropolis.ratio <- metropolis.ratio + propose.log.density(lambda.state, lambda.proposal)
    metropolis.ratio <- metropolis.ratio - propose.log.density(lambda.proposal, lambda.state)
    if ( log(runif(1)) < metropolis.ratio ) lambda.state <- lambda.proposal
  }
  lambda.state
}

sample.posterior$beta <- function ( data, param ) {
  prior.mean <- rep(0, 4)
  prior.variance <- diag(length(prior.mean)) * 30
  diagonal.prior.variance <- TRUE
  #'
  M <- diag(data$size) - param$lambda * param$W
  H <- 1:data$panel_length |> lapply( \ (time) cbind( data$X1[[time]], param$W %*% data$X2[[time]] ) )
  H.H <- H |> lapply( \ (x) t(x) %*% x ) |> Reduce(f=`+`)
  H.M.Y <- 1:data$panel_length |> lapply( \ (time) t(M %*% H[[time]]) %*% data$Y[[time]] ) |> Reduce(f=`+`)
  solved.prior.variance <- if ( diagonal.prior.variance ) {
    diag(length(prior.mean)) / prior.variance[1]
  } else {
    solve(prior.variance)
  }
  posterior.variance <- solve( solved.prior.variance + H.H / param$sigma2 )
  posterior.mean <- posterior.variance %*% ( solved.prior.variance %*% prior.mean + H.M.Y / param$sigma2 )
  rmvnorm <- function ( mean, cov ) as.vector( mean + t(chol(cov)) %*% rnorm(length(mean)) )
  rmvnorm(posterior.mean, posterior.variance)
}

sample.posterior$sigma2 <- function ( data, param ) {
  prior.a <- 1 # shape
  prior.b <- 1 # scale
  rinvgamma <- function ( shape, scale ) 1 / rgamma(1, shape=shape, scale=1/scale)
  posterior.a <- prior.a + 1/2 * data$size * data$panel_length
  posterior.b <- prior.b - LogLikelihoodKernelAll(
    data = data,
    lambda = param$lambda,
    W = param$W,
    beta1 = param$beta[1:2],
    beta2 = param$beta[3:4],
    sigma2 = 1 # NOTE: do not need sigma2
  )
  rinvgamma(posterior.a, posterior.b)
}

file.path("src", paste0("simulation-", TRIAL, "-net.cpp")) |> sourceCpp()
sample.posterior$net <- function ( data, param ) {
  SamplePosteriorNet(
    dimension = data$size^2,
    params = list(
      data = data,
      lambda = param$lambda,
      beta1 = param$beta[1:2],
      beta2 = param$beta[3:4],
      sigma2 = param$sigma2,
      theta = param$theta,
      size = data$size
    ),
    chain_length = 10,
    cycle_per_iteration = 6
  )
}

file.path("src", paste0("simulation-", TRIAL, "-theta.cpp")) |> sourceCpp()
sample.posterior$theta <- function ( data, param ) {
  SamplePosteriorTheta(
    param$net,
    parameters_dimension = 3,
    prior_mean = c(-3, 1, 1),
    prior_variance = 0.2,
    chain_length = 50
  ) |> as.vector()
}

### Gibbs #####################################################################

gibbs <- list()

gibbs$output.dir <- paste0("output-param-", TRIAL)
gibbs$save.name  <- function ( iteration ) file.path(gibbs$output.dir, paste0("param-", iteration, ".rds"))
gibbs$save.param <- function ( iteration, param ) saveRDS(param, file=gibbs$save.name(iteration))
gibbs$touch.output.dir <- function () {
  if ( !dir.exists(gibbs$output.dir) ) {
    dir.create(gibbs$output.dir)
  } else {
    paste0("Default output directory ", file.path(gibbs$output.dir), " already exists") |> message()
  }
}

gibbs$log.progress <- function ( iteration, param, chain.length ) {
  show.iteration <- sprintf( paste0("(Iteration %", as.integer(log10(chain.length))+1,"d)"), iteration )
  show.lambda    <- do.call( sprintf, list("lambda: %+1.5f", param$lambda) )
  show.sigma2    <- do.call( sprintf, list("sigma2: %2.2f",  param$sigma2) )
  show.beta      <- do.call( sprintf, c(list(paste0(c("beta:",  rep("%1.2f",length(param$beta))),  collapse=" ")), as.list(param$beta))  )
  show.theta     <- do.call( sprintf, c(list(paste0(c("theta:", rep("%1.2f",length(param$theta))), collapse=" ")), as.list(param$theta)) )
  cat(format(Sys.time(),"(%F %H:%M:%S)"), show.iteration, show.lambda, show.beta, show.sigma2, show.theta, fill = TRUE)
}

gibbs$run <- function (
  data,
  param = init.param(data),
  chain.length = 1e3,
  print.interval = 1,
  start.from.iteration = 1
) {
  gibbs$touch.output.dir()
  param <- init.param(data)
  for ( iteration in start.from.iteration:(chain.length-start.from.iteration+1) ) {
    param$beta   <- sample.posterior$beta(data,param)
    param$net    <- sample.posterior$net(data,param)
    param$W      <- RowNormalize(param$net)
    param$theta  <- sample.posterior$theta(data,param)
    param$lambda <- sample.posterior$lambda(data,param)
    param$sigma2 <- sample.posterior$sigma2(data,param)
    gibbs$save.param(iteration, param)
    if ( iteration %% print.interval == 0 ) gibbs$log.progress(iteration, param, chain.length)
  }
}

### Run Gibbs #################################################################

timer$tic()
gibbs$run(
  data = data,
  param = init.param(data),
  chain.length = 3e3,
  print.interval = 1,
  start.from.iteration = 1
)
timer$toc()
