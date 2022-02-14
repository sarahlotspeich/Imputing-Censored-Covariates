#' Generate data for linear regression with censored covariate
#'
#' Generates data.frame pertaining to the linear model, Y = beta0 + X betaX + Z betaZ + epsilon
#' and a censored observation of the covariate X by censoring variable C
#'
#' @param n sample size for single simulation
#' @param n.sims number of simulations
#' @param beta0 intercept
#' @param betaX coefficient for covariate \code{x}
#' @param betaZ coefficient for covariate(s) \code{z}. If \code{betaZ} NULL, then it is assumed that the model does not include \code{z}.
#' @param sigma sd for random error in outcome model
#' @param x censored covariate. If \code{x} NULL, then \code{x} is generated from Weibull(\code{xshape}, \code{xscale}).
#' @param z matrix of observed covariate(s). If \code{betaZ} is provided and \code{z} NULL, then \code{z} is generated from N_p(0, I).
#' @param xshape shape parameter for simulation of \code{x} (if \code{x} NULL)
#' @param xscale scale parameter for simulation of \code{x} (if \code{x} NULL)
#' @param c.lower minimum for simulation of \code{c}
#' @param c.upper maximum for simulation of \code{c}
#'
#' @return a data.frame with the following elements:
#' \itemize{
#' \item{subj.id} subject id
#' \item{sim.id} simulation id
#' \item{y} response
#' \item{x} true value of the covariate subject to censoring
#' \item{t} minimum of covariate \code{x} and censoring variable \code{c}
#' \item{event} 1 if \code{x} <= \code{c}; 0 otherwise
#' \item{z} observed covariates
#' }
#' 
#' @export
generate_data = function(n, n.sims = 1,
                         beta0, betaX, betaZ = NULL, sigma = 1,
                         x = NULL, z = NULL,
                         xshape = 0.75, xscale = 0.25,
                         c.lower = 0, c.upper = 1) {
  # test for bad input
  if (n %% 1 != 0 | n < 1) {stop("n must be positive integer")}
  if (n.sims %% 1 != 0 | n < 1) {stop("n.sims must be positive integer")}
  if (!is.null(betaZ) & !is.null(z)) {
    if (length(betaZ) != ncol(z)) {stop("number of parameters in betaZ must equal number of columns in z")}
  }
  if (xshape <= 0) {stop("xshape must be positive")}
  if (xscale <= 0) {stop("xscale must be positive")}
  if (c.lower >= c.upper) {stop("c.lower must be less than c.upper")}
    
  # total number of rows
  N <- n*n.sims

  # subject id and simulation id
  subj.id <- rep(1:n, n.sims)
  sim.id <- rep(1:n.sims, each = n)

  # Data generation
  # Generate x if not provided
  if (is.null(x)) { x <- rweibull(n = N, shape = xshape, scale = xscale) }

  # Generate random error and outcome
  epsilon <- rnorm(n = N, sd = sigma)
  y <- beta0 + betaX * x  + epsilon

  # Censoring mechanism
  cens <- runif(n = N, min = c.lower, max = c.upper)
  event <- as.numeric(x <= cens)
  t <- ifelse(event, x, cens)

  # If no z coefficient is provided, return id's, response, x, t, and event
  if (is.null(betaZ)) {
    data.frame(subj.id, sim.id, y, x, t, event) %>%
      return()
  }
  else {
    pZ <- length(betaZ)
    # Generate design matrix fof covariates z if not provided
    if (is.null(z)) { z <- rnorm(n = N*pZ, mean = 0, sd = 1) %>% matrix(nrow = N, ncol = pZ) }
    colnames(z) = paste0("Z", 1:pZ)
    # Add the effect of z to the outcome
    y <- y + z %*% betaZ
    # return id's, response, x, t, event, and z
    data.frame(subj.id, sim.id, y, x, t, event) %>%
      cbind(z) %>%
      return()
  }
}
