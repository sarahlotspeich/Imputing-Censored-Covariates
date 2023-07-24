mean_life = function(lin_pred, fit, dist) {
  if (dist %in% c("weibull", "exponential", "rayleigh")) {
    est_ml = exp(lin_pred) * gamma(1 + fit$scale)
  } #else if (dist %in% c("lognormal", "loggaussian")) {}
  return(est_ml)
}