#' Generate data for linear regression with censored covariate
#'
#' Generates data.frame pertaining to the linear model, Y = beta0 + X betaX + Z betaZ + epsilon
#' and a censored observation of the covariate X by censoring variable C
#'
#' @param n sample size for single simulation
#' @param n.sims number of simulations
#' @param beta0 intercept
#' @param betaX coefficient for covariate x
#' @param betaZ coefficient for covariate(s) z. If betaZ NULL, then it is assumed that the model does not include z.
#' @param sigma sd for random error in outcome model
#' @param x the censored covariate. If x NULL, then x is generated from Weibull(xshape, xscale).
#' @param z the observed covariate(s). If betaZ is provided and z NULL, then z is generated from N_p(0, I).
#' @param xshape shape parameter for simulation of x (if x NULL)
#' @param xscale scale parameter for simulation of x (if x NULL)
#' @param c.lower minimum for simulation of c
#' @param c.upper maximum for simulation of c
#'
#' @return a data.frame with the following elements:
#' \itemize{
#' \item{subj.id}
#' \item{sim.id}
#' \item{y}
#' \item{t}
#' \item{delta}
#' \item{z}
#' }
#'
#' @export
generate.data = function(n, n.sims = 1,
                         beta0, betaX, betaZ = NULL, sigma = 1,
                         x = NULL, z = NULL,
                         xshape = 0.75, xscale = 0.25,
                         c.lower = 0, c.upper = 1) {
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
  delta <- as.numeric(x <= cens)
  t <- ifelse(delta, x, cens)

  # If no z coefficient is provided, return id's, response, t, and delta
  if (is.null(betaZ)) {
    data.frame(subj.id, sim.id, y, t, delta) %>%
      return()
  }
  else {
    pZ <- length(betaZ)
    # Generate design matrix fof covariates z if not provided
    if (is.null(z)) { z <- rnorm(n = N*pZ, mean = 0, sd = 1) %>% matrix(nrow = N, ncol = pZ) }
    colnames(z) = paste0("Z", 1:pZ)
    # Add the effect of z to the outcome
    y <- y + z %*% betaZ
    # return id's, response, t, delta, and z
    data.frame(subj.id, sim.id, y, t, delta) %>%
      cbind(z) %>%
      return()
  }
}
