#' Generate data with a censored predictor
#'
#' Generate simulated data with a randomly right-censored predictor
#'
#' @param sample_size Number of independent observations to generate.
#' @param setting Choose from six simulation settings in the manuscript: \code{"PH_Light"}, \code{"PH_Heavy"}, \code{"NPH_Light"}, \code{"NPH_Heavy"}, \code{"Light"}, or \code{"Heavy"}.
#' @param cens_type Type of censoring (\code{"right"} or \code{"left"}). Default is \code{"right"}.
#' @param cens_lod (Optional) Limit of detection. Default is \code{cens_lod = NULL}, which corresponds to random censoring.
#'
#' @return A dataframe with \code{sample_size} rows, each containing one simulated observation comprised of the following columns:
#' \item{z}{A fully observed standard log-normal covariate.}
#' \item{x}{True values for the randomly right-censored predictor. This is fully observed and intended for empirical comparisons, though it would not be in practice.}
#' \item{t}{Observed values for the randomly right-censored predictor.}
#' \item{d}{Event indicator for \code{t} defined such that \code{d = 1} for uncensored subjects and \code{0} otherwise.}
#'
#' @export

generate_data <- function(sample_size, setting, cens_type = "right", cens_lod = NULL) {
  if (cens_type == "right") {
    if (is.null(cens_lod)) {
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
      } else if (setting == "Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25)
        c <- rweibull(n = sample_size, shape = 1, scale = 2)
        t <- pmin(x, c)
        d <- as.numeric(x <= c)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25)
        c <- rweibull(n = sample_size, shape = 1, scale = 0.35)
        t <- pmin(x, c)
        d <- as.numeric(x <= c)
        sdat <- data.frame(z, x, t, d)
      }
    } else {
      if (setting == "PH_Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05 * z)
        t <- pmin(x, cens_lod)
        d <- as.numeric(x <= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "PH_Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05  * z)
        t <- pmin(x, cens_lod)
        d <- as.numeric(x <= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "NPH_Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
        t <- pmin(x, cens_lod)
        d <- as.numeric(x <= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "NPH_Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
        t <- pmin(x, cens_lod)
        d <- as.numeric(x <= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting %in% c("Light", "Heavy")) { 
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25)
        t <- pmin(x, cens_lod)
        d <- as.numeric(x <= cens_lod)
        sdat <- data.frame(z, x, t, d)
      }
    }  
  } else if (cens_type == "left") {
    if (is.null(cens_lod)) {
      if (setting == "PH_Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05 * z)
        c <- rweibull(n = sample_size, shape = 1, scale = 2)
        t <- pmax(x, c)
        d <- as.numeric(x >= c)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "PH_Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05  * z)
        c <- rweibull(n = sample_size, shape = 1, scale = 0.35)
        t <- pmax(x, c)
        d <- as.numeric(x >= c)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "NPH_Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
        c <- rweibull(n = sample_size, shape = 1, scale = 2)
        t <- pmax(x, c)
        d <- as.numeric(x >= c)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "NPH_Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
        c <- rweibull(n = sample_size, shape = 1, scale = 0.35)
        t <- pmax(x, c)
        d <- as.numeric(x >= c)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25)
        c <- rweibull(n = sample_size, shape = 1, scale = 2)
        t <- pmax(x, c)
        d <- as.numeric(x >= c)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25)
        c <- rweibull(n = sample_size, shape = 1, scale = 0.35)
        t <- pmax(x, c)
        d <- as.numeric(x >= c)
        sdat <- data.frame(z, x, t, d)
      }
    } else {
      if (setting == "PH_Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05 * z)
        t <- pmax(x, cens_lod)
        d <- as.numeric(x >= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "PH_Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25 + 0.05  * z)
        t <- pmax(x, cens_lod)
        d <- as.numeric(x >= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "NPH_Light") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
        t <- pmax(x, cens_lod)
        d <- as.numeric(x >= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting == "NPH_Heavy") {
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6 + 0.9 * z, scale = 0.25)
        t <- pmax(x, cens_lod)
        d <- as.numeric(x >= cens_lod)
        sdat <- data.frame(z, x, t, d)
      } else if (setting %in% c("Light", "Heavy")) { 
        z <- rlnorm(n = sample_size)
        x <- rweibull(n = sample_size, shape = 0.6, scale = 0.25)
        t <- pmax(x, cens_lod)
        d <- as.numeric(x >= cens_lod)
        sdat <- data.frame(z, x, t, d)
      }
    }
  }
  return(sdat)
}