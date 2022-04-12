#' Custom conditional mean imputation (CMI) for a censored predictor with user-specified survival function
#'
#' Custom conditional mean imputation (CMI) for a censored predictor using (externally calculated) user-specified conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens).
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation.
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param useSURV Assumed survival function for \code{W} given \code{Z}. The only arguments to \code{useSURV} should be \code{W} and \code{Z}, in that order.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#'
#' @return
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_custom <- function(W, Delta, Z, data, useSURV, trapezoidal_rule = FALSE) {
  # Calculate survival with original model coefficients using custom function
  if(is.null(Z)) {
    data <- data.frame(data, surv = useSURV(data[, W]))
  } else {
    data <- data.frame(data, surv = useSURV(data[, W], data[, Z]))
  }

  # Order data by W
  data <- data[order(data[, W]), ]

  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1

  # Calculate imputed values
  data$imp <- data[, W]
  if (trapezoidal_rule) {
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, "surv")])

    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]

    # Censored subject values (to impute)
    t_cens <- data[data[, Delta] == 0, W]

    # Follow formula assuming AFT model for S(X|Z)
    for (x in which(!uncens)) {
      Cj <- data[x, W]
      Sj <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
      num <- sum((data_dist[-nrow(data_dist), W] >= Cj) * Sj * t_diff)
      denom <- data[x, "surv"]
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  } else {
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    if (is.null(Z)) {
      int_surv <- sapply(
        X = which(!uncens),
        FUN = function(i) {
          tryCatch(expr = integrate(f = function(t) useSURV(t), lower = data[i, W], upper = Inf)$value,
                   error = function(e) return(NA))
        }
      )
    } else {
      int_surv <- sapply(
        X = which(!uncens),
        FUN = function(i) {
          tryCatch(expr = integrate(f = function(t) useSURV(t, data[i, Z]), lower = data[i, W], upper = Inf)$value,
                   error = function(e) return(NA))
        }
      )
    }
    ## Calculate E(X|X>W,Z) = int_surv / surv(W|Z) + W
    data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }

  ## Check for infinite/NA imputed values
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] <- data[which(is.na(data$imp)), W]
  }
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }

  # Return input dataset with appended column imp containing imputed values
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
}
