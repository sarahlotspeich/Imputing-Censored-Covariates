#' Conditional mean imputation (CMI) for a censored predictor using the Kaplan-Meier estimator
#'
#' Nonparametric conditional mean imputation (CMI) for a censored predictor using a Kaplan-Meier product-limit estimator for estimate survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed binary covariate. If provided and \code{stratified = TRUE}, the Kaplan-Meier estimator will be stratified on \code{Z}.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param stratified A logical input for whether the Kaplan-Meier estimator should be stratified on \code{Z}. Default is \code{TRUE}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"drop-off"}, \code{"exponential"} (default), or \code{"weibull"}.
#'
#' @return A copy of \code{data} with added column \code{imp} containing the imputed values.
#'
#' @export
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom survival strata

cmi_km <- function(W, Delta, Z = NULL, data, stratified = TRUE, trapezoidal_rule = FALSE, surv_between = "carry-forward", surv_beyond = "exponential") {
  # Fit the Kaplan-Meier estimator for S(W)
  fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ 1"))
  
  # If (binary) Z is supplied, fit separate models for Z = 0/ Z = 1
  if (!is.null(Z) & stratified) {
    # Z = 0
    fit0 <- survfit(formula = fit_formula, 
                    data = data[data[, Z] == 0, ])
    surv_df0 <- data.frame(x = fit0$time, z = 0, surv = fit0$surv)
    colnames(surv_df0)[1:2] <- c(W, Z)
    
    # Z = 0
    fit1 <- survfit(formula = fit_formula, 
                    data = data[data[, Z] == 1, ])
    surv_df1 <- data.frame(x = fit1$time, z = 1, surv = fit1$surv)
    colnames(surv_df1)[1:2] <- c(W, Z)
    
    # Combine 
    surv_df <- rbind(surv_df0, surv_df1)
  } else {
    fit <- survfit(formula = fit_formula, 
                   data = data)
    surv_df <- data.frame(x = fit$time, surv = fit$surv)
    colnames(surv_df)[1] <- W
  }
  
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Save data frame with KM estimates (ordered)
  km_surv <- data[data[, Delta] == 1, ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1

  # Make survival at censored W NA so that we can interpolate/extrapolate
  data[!uncens, "surv"] <- NA
  
  # For people with events, obs = X
  data$imp <- data[, W]
  
  # If stratifying on Z, do the following within strata Z = 0 / Z = 1
  if (!is.null(Z) & stratified) {
    levelsZ <- unique(data[, Z])
    for (z in levelsZ) {
      # Assume survival at censored W < min(X) = 1
      minX <- data[min(which(uncens & data[, Z] == z)), W] 
      data[which(data[, W] < minX & data[, Z] == z), "surv"] <- 1
      
      # Interpolate baseline survival at censored W < \widetilde{X}
      Xtilde <- data[max(which(uncens & data[, Z] == z)), W] 
      needs_interp <- which(is.na(data[, "surv"]) & data[, Z] == z & data[, W] < Xtilde)
      data[needs_interp, "surv"] <- sapply(X = data[needs_interp, W], 
                                           FUN = interp_surv_between, 
                                           t = km_surv[which(km_surv[, Z] == z), W], 
                                           surv = km_surv[which(km_surv[, Z] == z), "surv"], 
                                           surv_between = surv_between)
      
      # Extrapolate baseline survival at censored W > \widetilde{X}
      if (surv_beyond == "weibull") {
        # Estimate Weibull parameters using constrained MLE
        Xtilde <- data[max(which(uncens & data[, Z] == z)), W] 
        SURVmax <- data[max(which(uncens & data[, Z] == z)), "surv"]
        weibull_params <- constr_weibull_mle(t = data[which(data[, Z] == z), W], 
                                             I_event = data[which(data[, Z] == z), Delta], 
                                             Xtilde = Xtilde, 
                                             rho = SURVmax, 
                                             alpha0 = 1E-4)
        
        # If weibull params don't converge, quit 
        if (any(is.na(weibull_params))) {
          data$imp <- NA
          return(list(imputed_data = data, code = FALSE))   
        }
        
        # If they do, extrapolate with them 
        needs_extrap <- which(!uncens & data[, W] > Xtilde & data[, Z] == z)
        data[needs_extrap, "surv"] <- sapply(X = data[needs_extrap, W], 
                                             FUN = extrap_surv_beyond, 
                                             t = km_surv[which(km_surv[, Z] == z), W], 
                                             surv = km_surv[which(km_surv[, Z] == z), "surv"], 
                                             surv_beyond = surv_beyond, 
                                             weibull_params = weibull_params)
      } else {
        needs_extrap <- which(!uncens & data[, W] > Xtilde & data[, Z] == z)
        data[needs_extrap, "surv"] <- sapply(X = data[needs_extrap, W], 
                                             FUN = extrap_surv_beyond, 
                                             t = km_surv[which(km_surv[, Z] == z), W], 
                                             surv = km_surv[which(km_surv[, Z] == z), "surv"], 
                                             surv_beyond = surv_beyond)
        weibull_params <- NULL
      }
      
      # Calculate imputed values E(X|X>W,Z)
      needs_impute <- which(!uncens & data[, Z] == z) # Row IDs in need of imputation with Z = z
      if (length(needs_impute) > 0) {
        if (trapezoidal_rule) {
          # Distinct rows (in case of non-unique obs values)
          data_dist <- unique(data[which(data[, Z] == z), c(W, Delta, Z, "surv")])
          
          # [T_{(i+1)} - T_{(i)}]
          t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
          
          # S(t+1) + S(t)
          surv_sum <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
          
          # Use trapezoidal approximation for integral
          for (i in needs_impute) {
            sum_surv_i <- sum((data_dist[- nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
            data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, "surv"])
          }
        } else {
          # Builds on the extend_surv function by raising S(t)
          to_integrate <- function(t) {
            surv <- sapply(X = t, 
                           FUN = extend_surv, 
                           t = km_surv[which(km_surv[, Z] == z), W], 
                           surv = km_surv[which(km_surv[, Z] == z), "surv"],
                           surv_between = surv_between, 
                           surv_beyond = surv_beyond, 
                           weibull_params = weibull_params)  
            surv
          }
          
          ## Use integrate() to approximate integral from W to \infty of S(t)
          int_surv <- sapply(
            X = needs_impute, 
            FUN = function(i) { 
              tryCatch(expr = integrate(f = to_integrate, 
                                        lower = data[i, W], 
                                        upper = Inf, 
                                        subdivisions = 2000)$value,
                       error = function(e) return(NA))
            }
          )
          
          ## Calculate E(X|X>W) = W + int_surv / surv(W)
          data$imp[needs_impute] <- data[needs_impute, W] + int_surv / data[needs_impute, "surv"]
        }  
      }
    }
  } else {
    # Assume survival at censored W < min(X) = 1
    minX <- data[min(which(uncens)), W] 
    data[which(data[, W] < minX), "surv"] <- 1
    
    # Interpolate baseline survival at censored W < \widetilde{X}
    Xtilde <- data[max(which(uncens)), W] 
    needs_interp <- which(is.na(data[, "surv"]) & data[, W] < Xtilde)
    data[needs_interp, "surv"] <- sapply(X = data[needs_interp, W], 
                                         FUN = interp_surv_between, 
                                         t = km_surv[, W], 
                                         surv = km_surv[, "surv"], 
                                         surv_between = surv_between)
    
    # Extrapolate baseline survival at censored W > \widetilde{X}
    if (surv_beyond == "weibull") {
      # Estimate Weibull parameters using constrained MLE
      SURVmax <- data[max(which(uncens)), "surv"]
      weibull_params <- constr_weibull_mle(t = data[, W], 
                                           I_event = data[, Delta], 
                                           Xtilde = Xtilde, 
                                           rho = SURVmax, 
                                           alpha0 = 1E-4)
      
      # If weibull params don't converge, quit 
      if (any(is.na(weibull_params))) {
        data$imp <- NA
        return(list(imputed_data = data, code = FALSE))   
      }
      
      # If they do, extrapolate with them 
      needs_extrap <- which(!uncens & data[, W] > Xtilde)
      data[needs_extrap, "surv"] <- sapply(X = data[needs_extrap, W], 
                                           FUN = extrap_surv_beyond, 
                                           t = km_surv[, W], 
                                           surv = km_surv[, "surv"], 
                                           surv_beyond = surv_beyond, 
                                           weibull_params = weibull_params)
    } else {
      needs_extrap <- which(!uncens & data[, W] > Xtilde)
      data[needs_extrap, "surv"] <- sapply(X = data[needs_extrap, W], 
                                           FUN = extrap_surv_beyond, 
                                           t = km_surv[, W], 
                                           surv = km_surv[, "surv"], 
                                           surv_beyond = surv_beyond)
      weibull_params <- NULL
    }
    
    # Calculate imputed values E(X|X>W)
    if (trapezoidal_rule) {
      # Distinct rows (in case of non-unique obs values)
      data_dist <- unique(data[, c(W, Delta, Z, "surv")])
      
      # [T_{(i+1)} - T_{(i)}]
      t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
      
      # S(t+1) + S(t)
      surv_sum <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
      
      # Use trapezoidal approximation for integral
      for (i in which(!uncens)) {
        sum_surv_i <- sum((data_dist[- nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
        data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, "surv"])
      }
    } else {
      # Builds on the extend_surv function by raising S(t)
      to_integrate <- function(t) {
        surv <- sapply(X = t, 
                       FUN = extend_surv, 
                       t = surv_df[, W], 
                       surv = surv_df[, "surv"], 
                       surv_between = surv_between, 
                       surv_beyond = surv_beyond, 
                       weibull_params = weibull_params)  
        surv
      }
      
      ## Use integrate() to approximate integral from W to \infty of S(t)
      int_surv <- sapply(
        X = which(!uncens), 
        FUN = function(i) { 
          tryCatch(expr = integrate(f = to_integrate, 
                                    lower = data[i, W], 
                                    upper = Inf, 
                                    subdivisions = 2000)$value,
                   error = function(e) return(NA))
        }
      )
      
      ## Calculate E(X|X>W) = W + int_surv / surv(W)
      data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
    }
  }
  
  ## Check for infinite/NA imputed values 
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] <- data[which(is.na(data$imp)), W]
  }
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(list(imputed_data = data, code = !any(is.na(data$imp))))
}
