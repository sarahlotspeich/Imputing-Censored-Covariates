#' Conditional mean imputation with bootstrap
#'
#' Imputes censored covariates with their conditional mean given censored value and additional covariates (where supplied).
#'
#' @param fit A \code{coxph} or \code{survfit} imputation model object.
#' @param obs String column name for the censored covariate.
#' @param delta String column name for the censoring indicator of the covariate.
#' @param addl_covar (Optional) string or vector of strings for the additional fully-observed covariates. Default is \code{NULL}.
#' @param data Datafrane containing columns \code{obs}, \code{delta}, and (if provided) \code{addl_covar}.
#' @param approx_beyond Choice of approximation used to extrapolate the survival function beyond the last observed covariate value. Default is \code{"expo"} for the exponential approximation. Other choices include \code{"zero"} or \code{"carryforward"}.
#' @param M number of bootstrap samples to be taken from \code{data}
#'
#' @return A list of \code{M} dataframes, each of which is sampled with replacement from \code{data} and then augmented with a column of imputed covariate values called \code{imp}.
#'
#' @importFrom dplyr filter
#' 
#' @export
condl_mean_impute_bootstrap <- function(fit, obs, delta, addl_covar = NULL, data, approx_beyond = "expo", M) {
  # list of imputed datasets to be returned
  imputed.datasets = list()
  n = nrow(data)
  
  # create M bootstrap samples of size n
  bs.data = data[sample(x = 1:(M*n), size = n, replace = T), ]
  bs.data$m = rep(1:M, each = n)
  
  # perform conditional mean imputation on each bootstrap sample
  for (m in 1:M) {
    imputed.datasets[m] <- dplyr::filter(bs.data, m == m) %>%
      condl_mean_impute(fit = fit, obs = obs, delta = delta, addl_covar = addl_covar,
                        data = data, approx_beyond = approx_beyond)
  }
  
  # return list of imputed datasets
  return(imputed.datasets)
}
