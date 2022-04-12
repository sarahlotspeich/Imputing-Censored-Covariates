#' Kernel-smoothed Nelson-Aalen estimator of conditional survival
#'
#' Uses the STRIDE implementation to estimate survival of \code{W} given \code{Z} nonparametrically using the kernel-smoothed Nelson-Aalen estimator. 
#'
#' @param t Value at which the survival function is to be estimated.
#' @param W Numeric vector of observed covariate values (including censored opens). 
#' @param Delta Numeric vector of censoring indicators to accompany \code{W}. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Numeric vector of additional fully observed (continuous) covariate to accompany \code{W} and \code{Delta}.
#'
#' @return Dataframe
#'
#' @export
#' @import stats

stride_surv = function(t, W, Delta, Z){
  h = bw.nrd(Z)
  bandwidth = h * length(Xi.use) ^ (- 1 / 10)
  K = Kern.FUN(Z, Z, bandwidth)
  tmpind = (W <= t) & (Delta == 1)
  tj = W[tmpind]
  kerni.1 = t(t(K))
  pihamyt0.tj.ss = sum.I(tj, "<=", W, Vi = t(kerni.1))
  dLamhat.tj.ss = t((kerni.1[, tmpind])) / pihamyt0.tj.ss
  ret = apply(dLamhat.tj.ss, 2, sum, na.rm = TRUE)
  S.return = exp(- ret)
  return(data.frame(W, Delta, Z, surv = S.return))
}