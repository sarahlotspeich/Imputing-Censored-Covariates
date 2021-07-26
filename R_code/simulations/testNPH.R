source("~sarahlotspeich/Dropbox/UNC/Sarah-Lotspeich/code/mi/impute_multivar.R")

sim_seed <- 918
set.seed(sim_seed)

PH <- FALSE
N <- 1000
beta <- 1
pZ <- 0.5

res <- data.frame()
for (scaleC in c(0.35)) {
    new_res <- data.frame(sim = 1:1000, scaleC, N, PH, 
        b0_atem_bres = NA, b1_atem_bres = NA, b2_atem_bres = NA, b0_se_atem_bres = NA, b1_se_atem_bres = NA, b2_se_atem_bres = NA,
        b0_atem_sf = NA, b1_atem_sf = NA, b2_atem_sf = NA, b0_se_atem_sf = NA, b1_se_atem_sf = NA, b2_se_atem_sf = NA, 
        b0_sarah_km = NA, b1_sarah_km = NA, b2_sarah_km = NA, b0_se_sarah_km = NA, b1_se_sarah_km = NA, b2_se_sarah_km = NA)
    for (s in 1:1000) {
        # Simulate data
        z <- rbinom(n = N, size = 1, prob = pZ)
        if (PH) {
          # Proportional hazard for X
          x <- rweibull(n = N, shape = 0.75, scale = 0.25)
        } else {
          # Non-proportional hazard for X
          x <- rweibull(n = N, shape = (0.6 + 0.9 * z), scale = 0.25)
        }
        e <- rnorm(n = N, mean = 0, sd = 1)
        y <- 1 + beta * x + 0.25 * z + e
        c <- rweibull(n = N, shape = 1, scale = scaleC) # for light or scale = 0.35 for heavy 
        d <- as.numeric(x <= c)
        t <- pmin(x, c)
        x[d == 0] <- NA
        sim_dat <- data.frame(y, x, t, z, d)
        # CMI - Atem et al. (2017)
        ## Using survfit() to estimate baseline survival
        cmi_sf <- cmi_multivar(data = sim_dat, outcome = "y", time = "t", event = "d", covar = "z", approx_last = "expo", est_S0 = "survfit", nonparam = FALSE, M = 20)
        new_res[s, c("b0_atem_sf", "b1_atem_sf", "b2_atem_sf")] <- cmi_sf$coeff
        new_res[s, c("b0_se_atem_sf", "b1_se_atem_sf", "b2_se_atem_sf")] <- cmi_sf$se
        ## Using hardcoded Breslow Estimator to estimate baseline survival
        cmi_bres <- cmi_multivar(data = sim_dat, outcome = "y", time = "t", event = "d", covar = "z", approx_last = "expo", est_S0 = "breslow", nonparam = FALSE, M = 20)
        new_res[s, c("b0_atem_bres", "b1_atem_bres", "b2_atem_bres")] <- cmi_bres$coeff
        new_res[s, c("b0_se_atem_bres", "b1_se_atem_bres", "b2_se_atem_bres")] <- cmi_bres$se
        # CMI(NP) - PROPOSED adaptation of Atem et al. (2017)
        cmi_km <- cmi_multivar(data = sim_dat, outcome = "y", time = "t", event = "d", covar = "z", approx_last = "expo", nonparam = TRUE, M = 20)  
        new_res[s, c("b0_sarah_km", "b1_sarah_km", "b2_sarah_km")] <- cmi_km$coeff
        new_res[s, c("b0_se_sarah_km", "b1_se_sarah_km", "b2_se_sarah_km")] <- cmi_km$se
    }
    res <- rbind(res, new_res)
    write.csv(res, paste0("testNPH_N", N, "_beta", ifelse(beta < 1, beta * 100, beta), "_pZ", pZ * 100, "_seed", sim_seed, ".csv"), row.names = F)
}
