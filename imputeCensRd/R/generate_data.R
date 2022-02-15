#' Generate data with a censored predictor
#'
#' Generate simulated data with a randomly right-censored predictor
#'
#' @param sample_size Number of independent observations to generate.
#' @param setting Choose from four simulation settings in the manuscript: \code{"PH_Light"}, \code{"PH_Heavy"}, \code{"NPH_Light"}, or \code{"NPH_Heavy"}.
#'
#' @return A dataframe with \code{sample_size} rows, each containing one simulated observation comprised of the following columns:
#' \item{z}{A fully observed standard log-normal covariate.}
#' \item{x}{True values for the randomly right-censored predictor. This is fully observed and intended for empirical comparisons, though it would not be in practice.}
#' \item{t}{Observed values for the randomly right-censored predictor.}
#' \item{d}{Event indicator for \code{t} defined such that \code{d = 1} for uncensored subjects and \code{0} otherwise.}
#'
#' @export

generate_data <- function(sample_size, setting) {
  if (setting == "PH_Light") {
    z <- rlnorm(n = sample_size)
    x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05 * z)
    c <- rweibull(n = sample_size, shape = 1, scale = 2)
    t <- pmin(x, c)
    d <- as.numeric(x <= c)
    sdat <- data.frame(z, x, t, d)
  } else if (setting == "PH_Heavy") {
    z <- rlnorm(n = sample_size)
    x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05  * z)
    c <- rweibull(n = sample_size, shape = 1, scale = 0.35)
    t <- pmin(x, c)
    d <- as.numeric(x <= c)
    sdat <- data.frame(z, x, t, d)
  } else if (setting == "NPH_Light") {
    z <- rlnorm(n = sample_size)
    x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
    c <- rweibull(n = sample_size, shape = 1, scale = 2)
    t <- pmin(x, c)
    d <- as.numeric(x <= c)
    sdat <- data.frame(z, x, t, d)
  } else if (setting == "NPH_Heavy") {
    z <- rlnorm(n = sample_size)
    x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
    c <- rweibull(n = sample_size, shape = 1, scale = 0.35)
    t <- pmin(x, c)
    d <- as.numeric(x <= c)
    sdat <- data.frame(z, x, t, d)
  }
  return(sdat)
}