#' Single, fully parametric conditional mean imputation for a censored covariate (Weibull distribution)
#'
#' Single, fully parametric conditional mean imputation for a censored covariate using an accelerated failure-time model with a Weibull distribution to estimate the conditional survival function.
#'
#' @param imputation_formula imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param W character, column name for observed values of the censored covariate 
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z character vector, column name(s) of additional fully observed covariate(s).
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param infinite_integral (optional) logical, if \code{infinite_integral = TRUE} (default) then conditional means are found by integrating from \code{W} to \code{Inf}, whereas if \code{infinite_integral = FALSE} they are found by subtracting the integral from \code{0} to \code{W} from the mean. If \code{infinite_integral = NA} instead, the analytical solutions are used to find the conditional means.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival survreg 
#' @importFrom survival Surv 
cmi_fp_weibull = function(imputation_formula, W, Delta, Z, data, infinite_integral = TRUE, maxiter = 100) {
  # Fit AFT imputation model for X ~ Z 
  fit = survreg(formula = imputation_formula, 
                data = data, 
                dist = "weibull", 
                maxiter = maxiter)
  
  # Initialize imputed values 
  data$imp = data[, W] ## start with imp = W
  
  # Calculate linear predictor for AFT imputation model
  #lp = fit$coefficients[1] +
  #  data.matrix(data[, Z]) %*% matrix(data = fit$coefficients[- 1], ncol = 1) ## linear predictors
  lp = fit$linear.predictors ## linear predictors
  
  # Transform parameters to agree with R's weibull parameterization 
  weib_shape = 1 / fit$scale
  weib_scale = exp(lp)
  
  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1
  
  if (is.na(infinite_integral)) {
    # Use closed-form to compute the conditional mean
    ## Transform parameters to agree with paper's parameterization
    alpha = weib_shape 
    lambda = weib_scale ^ (- weib_shape)
    
    # Use formula from Appendix for right-censored
    ## Save quantities for use in formula
    inside_exp = lambda[which(!uncens)] * data[which(!uncens), W] ^ alpha ## inside exp() for Weibull survival function
    gamma_surv = pgamma(q = inside_exp, 
                        shape = 1 / alpha, 
                        scale = 1, 
                        lower.tail = FALSE) ## survival function of a gamma
    data[which(!uncens), "imp"] = data[which(!uncens), W] * exp(- inside_exp) + 
      gamma(1 / alpha) / (alpha * lambda[which(!uncens)] ^ (1 / alpha)) * gamma_surv ## start with numerator
    data[which(!uncens), "imp"] = data[which(!uncens), "imp"] / exp(- inside_exp) ## divide by denominator
  } else {
    if (infinite_integral) {
      # Use adaptive quadrature to estimate
      ## integral from X = 0 to X = Wi
      int_surv_W_to_Inf = sapply(
        X = which(!uncens), 
        FUN = function(i) { 
          integrate(f = pweibull, 
                    lower = data[i, W], 
                    upper = Inf, 
                    shape = weib_shape, 
                    scale = weib_scale[i],
                    lower.tail = FALSE)$value
        }
      )
    } else if (!infinite_integral) {
      # Calculate mean life = integral from 0 to \infty of S(t|Z)
      est_ml = weib_scale[which(!uncens)] * gamma(1 + fit$scale) ## using formula for Weibull distribution
      
      # Use adaptive quadrature to estimate
      ## integral from X = 0 to X = Wi
      int_surv_0_to_W = sapply(
        X = which(!uncens), 
        FUN = function(i) { 
          integrate(f = pweibull, 
                    lower = 0, 
                    upper = data[i, W], 
                    shape = weib_shape, 
                    scale = weib_scale[i],
                    lower.tail = FALSE)$value
        }
      )
      
      # Subtract integral from mean life to get 
      ## integral from X = Wi to X = Infinity
      int_surv_W_to_Inf = est_ml - int_surv_0_to_W
    }
    
    ## Calculate S(W|Z)
    surv = pweibull(q = data[which(!uncens), W], 
                    shape = weib_shape, 
                    scale = weib_scale[which(!uncens)], 
                    lower.tail = FALSE)
    
    ## Calculate MRL(W) = int_surv / S(W|Z)
    est_mrl = int_surv_W_to_Inf / surv
    
    ## Calculate imputed value E(X|X>W,Z) = W + MRL(W)
    data[which(!uncens), "imp"] = data[which(!uncens), W] + est_mrl
    
    ## Check for infinite/NA imputed values 
    # if (any(is.na(data$imp))) {
    #   data$imp[which(is.na(data$imp))] = data[which(is.na(data$imp)), W]
    # }
    # if (any(data$imp == Inf)) {
    #   data$imp[which(data$imp == Inf)] = data[which(data$imp == Inf), W]
    # }
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(list(imputed_data = data, code = !any(is.na(data$imp))))  
}