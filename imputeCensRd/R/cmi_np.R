#' Nonparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Nonparametric conditional mean imputation (CMI) for a censored predictor using a modified Nelson-Aalen estimator (STRIDE) to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param approx_beyond Choice of approximation used to extrapolate the survival function beyond the last observed covariate value. Default is \code{"expo"} for the exponential approximation. Other choices include \code{"zero"} or \code{"carryforward"}.
#'
#' @return A copy of \code{data} with added column \code{imp} containing the imputed values.
#'
#' @export

cmi_np <- function(W, Delta, Z, data, approx_beyond = "expo") {
  # Estimate survival from kernel-smoothed Nelson-Aalen estimator (STRIDE)
  ## Create artificial subgroup/categorical covariates 
  data$subgroup <- data$cat <- 1
  ## Reorder columns per STRIDE documentation
  data <- data[, c(W, Delta, "cat", Z, "subgroup")]
  all_t <- data[, "t"]
  surv_df <- data[, c("cat", Z, "subgroup")]
  surv_df$surv <- NA
  for (ind in 1:length(all_t)) {
    if (ind == length(all_t)) {
      stride_fit <- S.NPNA.zw(t = all_t[ind], data = data, newdata = surv_df[(ind - 1):ind, ], weight = NULL)
      surv_df$surv[ind] <- stride_fit$data.out$survival.v.new[2]
    } else {
      stride_fit <- S.NPNA.zw(t = all_t[ind], data = data, newdata = surv_df[ind:(ind + 1), ], weight = NULL)
      surv_df$surv[ind] <- stride_fit$data.out$survival.v.new[1]
    }
  }
  surv_df$t <- all_t
  
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  # Order data by W
  data <- data[order(data[, W]), ]
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  if (any(is.na(data[, "surv"]))) {
    # For censored subjects, survival is average of times right before/after
    suppressWarnings(
      data[is.na(data[, "surv"]), "surv"] <- sapply(X = data[is.na(data[, "surv"]), W], FUN = impute_censored_surv, time = W, event = Delta, surv = "surv", data = data)
    )
  }
  
  # Assume survival = 1 before first observed covariate value
  if (any(data[, W] < min(data[uncens, W]))) {
    cens_before <- which(data[, W] < min(data[uncens, W]))
    data[cens_before, "surv"] <- 1
  }
  
  # Extrapolate survival beyond last observed covariate
  if (any(data[, W] > max(data[uncens, W]))) {
    cens_after <- which(data[, W] > max(data[uncens, W]))
    t_cens_after <- data[cens_after, W]
    last_event_surv <- data[max(which(uncens)), "surv"]
    # Efron (1967) immediately goes to zero
    if (approx_beyond == "zero") { data[cens_after, "surv"] <- 0 }
    # Gill (1980) carry forward survival at last Delta
    if (approx_beyond == "carryforward") { data[cens_after, "surv"] <- last_event_surv }
    # Brown, Hollander, and Kowar (1974) exponential approx
    if (approx_beyond == "expo") {
      max_event <- max(which(data[, Delta] == 1))
      t_max_event <- data[max_event, W]
      surv_max_event <- data[max_event, "surv"]
      data[cens_after, "surv"] <- exp(t_cens_after * log(surv_max_event) / t_max_event)
    }
  }
  # Distinct rows (in case of non-unique W values)
  data_dist <- unique(data[, c(W, Delta, Z, "surv")])
  # [T_{(i+1)} - T_{(i)}]
  t_diff <- data_dist[-1, W] - data_dist[-nrow(data_dist), W]
  # Censored subject values (to impute)
  t_cens <- data[data[, Delta] == 0, W]
  # For people with events, W = X
  data$imp <- data[, W]
  
  # Follow formula without additional covariates Z
  for (x in which(!uncens)) {
    Cj <- data[x, W]
    Sj <- data_dist[-1, "surv"] + data_dist[-nrow(data_dist), "surv"]
    num <- sum((data_dist[- nrow(data_dist), W] >= Cj) * Sj * t_diff)
    denom <- data[x, "surv"]
    data$imp[x] <- (1 / 2) * (num / denom) + Cj
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(data)
}