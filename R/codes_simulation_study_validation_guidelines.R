

###############################
#### GUIDELINES VALIDATION ####
###############################

#### DO NOT RUN: COMPUTATIONALLY VERY EXPENSIVE ####

library(doParallel)
library(foreach)
library(progressr)
library(future.apply)

##############################
#### GRID = F/T | N = 100 ####
##############################

curve_data_gridF_n100
median_comb_gridF_n100 <- data.frame(mult = curve_data_gridF_n100$mult, 
                                     ncube = round(curve_data_gridF_n100$median),
                                     grid = rep(FALSE, nrow(curve_data_gridF_n100)))

curve_data_gridT_n100
median_comb_gridT_n100 <- data.frame(mult = curve_data_gridT_n100$mult, 
                                     ncube = round(curve_data_gridT_n100$median),
                                     grid = rep(TRUE, nrow(curve_data_gridT_n100)))

median_comb_n100 <- rbind(median_comb_gridF_n100, median_comb_gridT_n100)

list_all_scen_n100 <- list(
  inh_scen1_n100, inh_scen2_n100, inh_scen3_n100, inh_scen4_n100,
  inh_scen5_n100, inh_scen6_n100, inh_scen7_n100, inh_scen8_n100
)

list_ppp_all_scenario_gridF_n100 <- list()

seed <- 2
for(i in seq_along(list_all_scen_n100)) {
  set.seed(seed + i)
  id <- sample(1:100, 10) 
  list_ppp_all_scenario_gridF_n100 <- c(
    list_ppp_all_scenario_gridF_n100,
    list_all_scen_n100[[i]][id]
  )
}

id_scenario <- rep(1:8, each = 10)

n_scenarios <- 8
n_processi_per_scenario <- 10
n_comb <- nrow(median_comb_n100)

dataset_n100_pvalue <- do.call(rbind, lapply(1:n_scenarios, function(s) {
  do.call(rbind, lapply(1:n_processi_per_scenario, function(p) {
    df <- median_comb_n100
    df$id_scenario <- s
    df$id_processo <- p
    df$pvalue <- NA_real_
    df[, c("id_scenario", "id_processo", "mult", "ncube", "grid", "pvalue")]
  }))
}))

n_comb_total <- nrow(dataset_n100_pvalue)

plan(sequential)
plan(multisession, workers = 60)
handlers(global = TRUE)

calc_pvalue <- function(row_index, p) {
  row <- dataset_n100_pvalue[row_index, ]
  idx <- (row$id_scenario - 1L) * n_processi_per_scenario + row$id_processo
  ppp_true <- list_ppp_all_scenario_gridF_n100[[idx]]
  current_scenario <- row$id_scenario
  comb <- row
  
  result <- tryCatch({
    switch(as.character(current_scenario),
           
           "1" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x,
                                      mult = comb$mult, ncube = comb$ncube, grid = comb$grid, seed = 2)
             est_par <- c(as.numeric(est$IntCoefs[1]), as.numeric(est$IntCoefs[2]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp(lambda = function(x, y, t, a) exp(a[1] + a[2]*x),
                              par = est_par, nsim = 50)
             Khat_sim <- matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp))
             for (j in seq_along(sim_ppp)) {
               pred_int <- exp(est_par[1] + est_par[2]*sim_ppp[[j]]$df$x)
               Khat_sim_j <- Khat_spatial.3D(sim_ppp[[j]], lambda = pred_int,
                                             correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat,
                                                sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "2" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x  + a[3] * y  + a[4] * t)},
                              par = c(est_par[1], est_par[2], est_par[3], est_par[4]), nsim = 50)
             
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$df$x +
                                 est_par[3] * sim_ppp[[j]]$df$y + est_par[4] * sim_ppp[[j]]$df$t)
               Khat_sim_j <- Khat_spatial.3D(sim_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "3" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ cov1, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                                 par = c(est_par[1], est_par[2]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$V1)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "4" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ cov1 + cov2 + cov3, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, verbose = F,
                                                 par = c(est_par[1], est_par[2], est_par[3], est_par[4]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$V1 +
                                 est_par[3] * sim_ppp[[j]]$V2 + est_par[4] * sim_ppp[[j]]$V3)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "5" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + cov1 + cov2 + cov3, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]),
                          est_param5 = as.numeric(est$IntCoefs[5]), est_param6 = as.numeric(est$IntCoefs[6]),
                          est_param7 = as.numeric(est$IntCoefs[7]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                                 par = c(est_par[1], est_par[2], est_par[3], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$x + est_par[3] * sim_ppp[[j]]$y + est_par[4] * sim_ppp[[j]]$t +
                                 est_par[5] * sim_ppp[[j]]$V1 + est_par[6] * sim_ppp[[j]]$V2 + est_par[7] * sim_ppp[[j]]$V3)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "6" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[2]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[3]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[1], nsim = 50)
             set.seed(2)
             g2_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[2], nsim = 50)
             set.seed(3)
             g3_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[3], nsim = 50)
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- stpm(data.frame(x = c(g1_sim_ppp[[si]]$df$x, g2_sim_ppp[[si]]$df$x, g3_sim_ppp[[si]]$df$x),
                                                y = c(g1_sim_ppp[[si]]$df$y, g2_sim_ppp[[si]]$df$y, g3_sim_ppp[[si]]$df$y),
                                                t = c(g1_sim_ppp[[si]]$df$t, g2_sim_ppp[[si]]$df$t, g3_sim_ppp[[si]]$df$t),
                                                m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]]$df)), rep("b", nrow(g2_sim_ppp[[si]]$df)), rep("c", nrow(g3_sim_ppp[[si]]$df))))))
             }
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$df$m1
               linear_predictor_sim <- param_vec[marks_vec_sim]
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "7" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[5]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
                          est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), est_param6 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[1], est_par[4], est_par[5], est_par[6]), nsim = 50)
             set.seed(2)
             g2_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[2], est_par[4], est_par[5], est_par[6]), nsim = 50)
             set.seed(3)
             g3_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[3], est_par[4], est_par[5], est_par[6]), nsim = 50)
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- stpm(data.frame(x = c(g1_sim_ppp[[si]]$df$x, g2_sim_ppp[[si]]$df$x, g3_sim_ppp[[si]]$df$x),
                                                y = c(g1_sim_ppp[[si]]$df$y, g2_sim_ppp[[si]]$df$y, g3_sim_ppp[[si]]$df$y),
                                                t = c(g1_sim_ppp[[si]]$df$t, g2_sim_ppp[[si]]$df$t, g3_sim_ppp[[si]]$df$t),
                                                m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]]$df)), rep("b", nrow(g2_sim_ppp[[si]]$df)), rep("c", nrow(g3_sim_ppp[[si]]$df))))))
             }
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$df$m1
               linear_predictor_sim <- param_vec[marks_vec_sim] + est_par[4] * sim_ppp[[j]]$df$x + est_par[5] * sim_ppp[[j]]$df$y + est_par[6] * sim_ppp[[j]]$df$t
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "8" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[8]),
                          est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), 
                          est_param6 = as.numeric(est$IntCoefs[4]), est_param7 = as.numeric(est$IntCoefs[5]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[1], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             set.seed(2)
             g2_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[2], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             set.seed(3)
             g3_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[3], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- data.frame(x = c(g1_sim_ppp[[si]]$x, g2_sim_ppp[[si]]$x, g3_sim_ppp[[si]]$x),
                                           y = c(g1_sim_ppp[[si]]$y, g2_sim_ppp[[si]]$y, g3_sim_ppp[[si]]$y),
                                           t = c(g1_sim_ppp[[si]]$t, g2_sim_ppp[[si]]$t, g3_sim_ppp[[si]]$t),
                                           cov1 = c(g1_sim_ppp[[si]]$V1, g2_sim_ppp[[si]]$V1, g3_sim_ppp[[si]]$V1),
                                           m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]])), rep("b", nrow(g2_sim_ppp[[si]])), rep("c", nrow(g3_sim_ppp[[si]])))))
             }
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$m1
               linear_predictor_sim <- param_vec[marks_vec_sim] + est_par[4] * sim_ppp[[j]]$x + est_par[5] * sim_ppp[[j]]$y + est_par[6] * sim_ppp[[j]]$t + est_par[7] * sim_ppp[[j]]$cov1
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           }
    )
  }, error = function(e) NA_real_)
  
  p() 
  return(result)
}

with_progress({
  p <- progressor(steps = n_comb_total)
  results_n100_pvalue <- future_lapply(1:n_comb_total, calc_pvalue, p = p)
})

dataset_n100_pvalue$pvalue <- unlist(results_n100_pvalue)
save(dataset_n100_pvalue, file = "dataset_n100_pvalue_parallel.RData")


dataset_n100_pvalue_scen678 <- dataset_n100_pvalue %>% filter(id_scenario %in% c(6,7,8))
summary(dataset_n100_pvalue_scen678)

list_ppp_all_scenario_gridF_n100_scen678 <- list_ppp_all_scenario_gridF_n100[51:80]

n_comb_scen678 <- nrow(dataset_n100_pvalue_scen678)

with_progress({
  p <- progressor(steps = n_comb_scen678)
  results_n100_pvalue_scen678 <- future_lapply(1:n_comb_scen678, calc_pvalue, p = p)
})


dataset_n100_pvalue_scen678$pvalue <- unlist(results_n100_pvalue_scen678)
summary(dataset_n100_pvalue_scen678)

dataset_n100_pvalue_final <- dataset_n100_pvalue

idx_replace <- dataset_n100_pvalue_final$id_scenario %in% c(6, 7, 8)

dataset_n100_pvalue_final <- dataset_n100_pvalue_final[order(dataset_n100_pvalue_final$id_scenario,
                                                             dataset_n100_pvalue_final$id_processo,
                                                             dataset_n100_pvalue_final$mult,
                                                             dataset_n100_pvalue_final$grid), ]

dataset_n100_pvalue_scen678 <- dataset_n100_pvalue_scen678[order(dataset_n100_pvalue_scen678$id_scenario,
                                                                 dataset_n100_pvalue_scen678$id_processo,
                                                                 dataset_n100_pvalue_scen678$mult,
                                                                 dataset_n100_pvalue_scen678$grid), ]

dataset_n100_pvalue_final[idx_replace, ] <- dataset_n100_pvalue_scen678
length(which(dataset_n100_pvalue_final$pvalue < 0.05))
save(dataset_n100_pvalue_final, file = "dataset_n100_pvalue_final.RData")





##############################
#### GRID = F/T | N = 250 ####
##############################

curve_data_gridF_n250
median_comb_gridF_n250 <- data.frame(mult = curve_data_gridF_n250$mult, 
                                     ncube = round(curve_data_gridF_n250$median),
                                     grid = rep(FALSE, nrow(curve_data_gridF_n250)))

curve_data_gridT_n250
median_comb_gridT_n250 <- data.frame(mult = curve_data_gridT_n250$mult, 
                                     ncube = round(curve_data_gridT_n250$median),
                                     grid = rep(TRUE, nrow(curve_data_gridT_n250)))

median_comb_n250 <- rbind(median_comb_gridF_n250, median_comb_gridT_n250)

list_all_scen_n250 <- list(
  inh_scen1_n250, inh_scen2_n250, inh_scen3_n250, inh_scen4_n250,
  inh_scen5_n250, inh_scen6_n250, inh_scen7_n250, inh_scen8_n250
)

list_ppp_all_scenario_gridF_n250 <- list()

seed <- 2
for(i in seq_along(list_all_scen_n250)) {
  set.seed(seed + i)
  id <- sample(1:100, 10) 
  list_ppp_all_scenario_gridF_n250 <- c(
    list_ppp_all_scenario_gridF_n250,
    list_all_scen_n250[[i]][id]
  )
}

id_scenario <- rep(1:8, each = 10)

n_scenarios <- 8
n_processi_per_scenario <- 10
n_comb <- nrow(median_comb_n250)

dataset_n250_pvalue <- do.call(rbind, lapply(1:n_scenarios, function(s) {
  do.call(rbind, lapply(1:n_processi_per_scenario, function(p) {
    df <- median_comb_n250
    df$id_scenario <- s
    df$id_processo <- p
    df$pvalue <- NA_real_
    df[, c("id_scenario", "id_processo", "mult", "ncube", "grid", "pvalue")]
  }))
}))

n_comb_total <- nrow(dataset_n250_pvalue)

plan(sequential)
plan(multisession, workers = 60)
handlers(global = TRUE)

calc_pvalue_n250 <- function(row_index, p) {
  row <- dataset_n250_pvalue[row_index, ]
  idx <- (row$id_scenario - 1L) * n_processi_per_scenario + row$id_processo
  ppp_true <- list_ppp_all_scenario_gridF_n250[[idx]]
  current_scenario <- row$id_scenario
  comb <- row
  
  result <- tryCatch({
    switch(as.character(current_scenario),
           
           "1" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x,
                                      mult = comb$mult, ncube = comb$ncube, grid = comb$grid, seed = 2)
             est_par <- c(as.numeric(est$IntCoefs[1]), as.numeric(est$IntCoefs[2]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp(lambda = function(x, y, t, a) exp(a[1] + a[2]*x),
                              par = est_par, nsim = 50)
             Khat_sim <- matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp))
             for (j in seq_along(sim_ppp)) {
               pred_int <- exp(est_par[1] + est_par[2]*sim_ppp[[j]]$df$x)
               Khat_sim_j <- Khat_spatial.3D(sim_ppp[[j]], lambda = pred_int,
                                             correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat,
                                                sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "2" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x  + a[3] * y  + a[4] * t)},
                              par = c(est_par[1], est_par[2], est_par[3], est_par[4]), nsim = 50)
             
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$df$x +
                                 est_par[3] * sim_ppp[[j]]$df$y + est_par[4] * sim_ppp[[j]]$df$t)
               Khat_sim_j <- Khat_spatial.3D(sim_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "3" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ cov1, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                                 par = c(est_par[1], est_par[2]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$V1)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "4" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ cov1 + cov2 + cov3, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, verbose = F,
                                                 par = c(est_par[1], est_par[2], est_par[3], est_par[4]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$V1 +
                                 est_par[3] * sim_ppp[[j]]$V2 + est_par[4] * sim_ppp[[j]]$V3)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "5" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + cov1 + cov2 + cov3, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]),
                          est_param5 = as.numeric(est$IntCoefs[5]), est_param6 = as.numeric(est$IntCoefs[6]),
                          est_param7 = as.numeric(est$IntCoefs[7]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                                 par = c(est_par[1], est_par[2], est_par[3], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$x + est_par[3] * sim_ppp[[j]]$y + est_par[4] * sim_ppp[[j]]$t +
                                 est_par[5] * sim_ppp[[j]]$V1 + est_par[6] * sim_ppp[[j]]$V2 + est_par[7] * sim_ppp[[j]]$V3)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "6" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[2]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[3]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[1], nsim = 50)
             set.seed(2)
             g2_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[2], nsim = 50)
             set.seed(3)
             g3_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[3], nsim = 50)
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- stpm(data.frame(x = c(g1_sim_ppp[[si]]$df$x, g2_sim_ppp[[si]]$df$x, g3_sim_ppp[[si]]$df$x),
                                                y = c(g1_sim_ppp[[si]]$df$y, g2_sim_ppp[[si]]$df$y, g3_sim_ppp[[si]]$df$y),
                                                t = c(g1_sim_ppp[[si]]$df$t, g2_sim_ppp[[si]]$df$t, g3_sim_ppp[[si]]$df$t),
                                                m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]]$df)), rep("b", nrow(g2_sim_ppp[[si]]$df)), rep("c", nrow(g3_sim_ppp[[si]]$df))))))
             }
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$df$m1
               linear_predictor_sim <- param_vec[marks_vec_sim]
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "7" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[5]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
                          est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), est_param6 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[1], est_par[4], est_par[5], est_par[6]), nsim = 50)
             set.seed(2)
             g2_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[2], est_par[4], est_par[5], est_par[6]), nsim = 50)
             set.seed(3)
             g3_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[3], est_par[4], est_par[5], est_par[6]), nsim = 50)
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- stpm(data.frame(x = c(g1_sim_ppp[[si]]$df$x, g2_sim_ppp[[si]]$df$x, g3_sim_ppp[[si]]$df$x),
                                                y = c(g1_sim_ppp[[si]]$df$y, g2_sim_ppp[[si]]$df$y, g3_sim_ppp[[si]]$df$y),
                                                t = c(g1_sim_ppp[[si]]$df$t, g2_sim_ppp[[si]]$df$t, g3_sim_ppp[[si]]$df$t),
                                                m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]]$df)), rep("b", nrow(g2_sim_ppp[[si]]$df)), rep("c", nrow(g3_sim_ppp[[si]]$df))))))
             }
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$df$m1
               linear_predictor_sim <- param_vec[marks_vec_sim] + est_par[4] * sim_ppp[[j]]$df$x + est_par[5] * sim_ppp[[j]]$df$y + est_par[6] * sim_ppp[[j]]$df$t
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "8" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[8]),
                          est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), 
                          est_param6 = as.numeric(est$IntCoefs[4]), est_param7 = as.numeric(est$IntCoefs[5]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[1], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             set.seed(2)
             g2_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[2], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             set.seed(3)
             g3_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[3], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- data.frame(x = c(g1_sim_ppp[[si]]$x, g2_sim_ppp[[si]]$x, g3_sim_ppp[[si]]$x),
                                           y = c(g1_sim_ppp[[si]]$y, g2_sim_ppp[[si]]$y, g3_sim_ppp[[si]]$y),
                                           t = c(g1_sim_ppp[[si]]$t, g2_sim_ppp[[si]]$t, g3_sim_ppp[[si]]$t),
                                           cov1 = c(g1_sim_ppp[[si]]$V1, g2_sim_ppp[[si]]$V1, g3_sim_ppp[[si]]$V1),
                                           m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]])), rep("b", nrow(g2_sim_ppp[[si]])), rep("c", nrow(g3_sim_ppp[[si]])))))
             }
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$m1
               linear_predictor_sim <- param_vec[marks_vec_sim] + est_par[4] * sim_ppp[[j]]$x + est_par[5] * sim_ppp[[j]]$y + est_par[6] * sim_ppp[[j]]$t + est_par[7] * sim_ppp[[j]]$cov1
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           }
    )
  }, error = function(e) NA_real_)
  
  p() 
  return(result)
}

with_progress({
  p <- progressor(steps = n_comb_total)
  results_n250_pvalue <- future_lapply(1:n_comb_total, calc_pvalue_n250, p = p)
})
dataset_n250_pvalue$pvalue <- unlist(results_n250_pvalue)
save(dataset_n250_pvalue, file = "dataset_n250_pvalue_parallel.RData")

dataset_n250_pvalue_na <- dataset_n250_pvalue[which(is.na(dataset_n250_pvalue$pvalue)),]
summary(dataset_n250_pvalue_na)
table(dataset_n250_pvalue_na$id_scenario)



##############################
#### GRID = F/T | N = 500 ####
##############################

curve_data_gridF_n500
median_comb_gridF_n500 <- data.frame(mult = curve_data_gridF_n500$mult, 
                                     ncube = round(curve_data_gridF_n500$median),
                                     grid = rep(FALSE, nrow(curve_data_gridF_n500)))

curve_data_gridT_n500
median_comb_gridT_n500 <- data.frame(mult = curve_data_gridT_n500$mult, 
                                     ncube = round(curve_data_gridT_n500$median),
                                     grid = rep(TRUE, nrow(curve_data_gridT_n500)))

median_comb_n500 <- rbind(median_comb_gridF_n500, median_comb_gridT_n500)

list_all_scen_n500 <- list(
  inh_scen1_n500, inh_scen2_n500, inh_scen3_n500, inh_scen4_n500,
  inh_scen5_n500, inh_scen6_n500, inh_scen7_n500, inh_scen8_n500
)

list_ppp_all_scenario_gridF_n500 <- list()

seed <- 2
for(i in seq_along(list_all_scen_n500)) {
  set.seed(seed + i)
  id <- sample(1:100, 10) 
  list_ppp_all_scenario_gridF_n500 <- c(
    list_ppp_all_scenario_gridF_n500,
    list_all_scen_n500[[i]][id]
  )
}

id_scenario <- rep(1:8, each = 10)

n_scenarios <- 8
n_processi_per_scenario <- 10
n_comb <- nrow(median_comb_n500)

dataset_n500_pvalue <- do.call(rbind, lapply(1:n_scenarios, function(s) {
  do.call(rbind, lapply(1:n_processi_per_scenario, function(p) {
    df <- median_comb_n500
    df$id_scenario <- s
    df$id_processo <- p
    df$pvalue <- NA_real_
    df[, c("id_scenario", "id_processo", "mult", "ncube", "grid", "pvalue")]
  }))
}))

n_comb_total <- nrow(dataset_n500_pvalue)

plan(sequential)
plan(multisession, workers = 60)
handlers(global = TRUE)

calc_pvalue_n500 <- function(row_index, p) {
  row <- dataset_n500_pvalue[row_index, ]
  idx <- (row$id_scenario - 1L) * n_processi_per_scenario + row$id_processo
  ppp_true <- list_ppp_all_scenario_gridF_n500[[idx]]
  current_scenario <- row$id_scenario
  comb <- row
  
  result <- tryCatch({
    switch(as.character(current_scenario),
           
           "1" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x,
                                      mult = comb$mult, ncube = comb$ncube, grid = comb$grid, seed = 2)
             est_par <- c(as.numeric(est$IntCoefs[1]), as.numeric(est$IntCoefs[2]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp(lambda = function(x, y, t, a) exp(a[1] + a[2]*x),
                              par = est_par, nsim = 50)
             Khat_sim <- matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp))
             for (j in seq_along(sim_ppp)) {
               pred_int <- exp(est_par[1] + est_par[2]*sim_ppp[[j]]$df$x)
               Khat_sim_j <- Khat_spatial.3D(sim_ppp[[j]], lambda = pred_int,
                                             correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat,
                                                sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "2" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x  + a[3] * y  + a[4] * t)},
                              par = c(est_par[1], est_par[2], est_par[3], est_par[4]), nsim = 50)
             
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$df$x +
                                 est_par[3] * sim_ppp[[j]]$df$y + est_par[4] * sim_ppp[[j]]$df$t)
               Khat_sim_j <- Khat_spatial.3D(sim_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "3" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ cov1, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                                 par = c(est_par[1], est_par[2]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$V1)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "4" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ cov1 + cov2 + cov3, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, verbose = F,
                                                 par = c(est_par[1], est_par[2], est_par[3], est_par[4]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$V1 +
                                 est_par[3] * sim_ppp[[j]]$V2 + est_par[4] * sim_ppp[[j]]$V3)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "5" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + cov1 + cov2 + cov3, mult = comb$mult, ncube = comb$ncube, grid = comb$grid, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
                          est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]),
                          est_param5 = as.numeric(est$IntCoefs[5]), est_param6 = as.numeric(est$IntCoefs[6]),
                          est_param7 = as.numeric(est$IntCoefs[7]))
             Khat_data <- Khat_spatial.3D(ppp_true, lambda = est$l, correction = "translate")
             sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                                 par = c(est_par[1], est_par[2], est_par[3], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               pred_int <- exp(est_par[1] + est_par[2] * sim_ppp[[j]]$x + est_par[3] * sim_ppp[[j]]$y + est_par[4] * sim_ppp[[j]]$t +
                                 est_par[5] * sim_ppp[[j]]$V1 + est_par[6] * sim_ppp[[j]]$V2 + est_par[7] * sim_ppp[[j]]$V3)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "6" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[2]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[3]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[1], nsim = 50)
             set.seed(2)
             g2_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[2], nsim = 50)
             set.seed(3)
             g3_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = est_par[3], nsim = 50)
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- stpm(data.frame(x = c(g1_sim_ppp[[si]]$df$x, g2_sim_ppp[[si]]$df$x, g3_sim_ppp[[si]]$df$x),
                                                y = c(g1_sim_ppp[[si]]$df$y, g2_sim_ppp[[si]]$df$y, g3_sim_ppp[[si]]$df$y),
                                                t = c(g1_sim_ppp[[si]]$df$t, g2_sim_ppp[[si]]$df$t, g3_sim_ppp[[si]]$df$t),
                                                m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]]$df)), rep("b", nrow(g2_sim_ppp[[si]]$df)), rep("c", nrow(g3_sim_ppp[[si]]$df))))))
             }
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$df$m1
               linear_predictor_sim <- param_vec[marks_vec_sim]
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "7" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[5]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
                          est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), est_param6 = as.numeric(est$IntCoefs[4]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[1], est_par[4], est_par[5], est_par[6]), nsim = 50)
             set.seed(2)
             g2_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[2], est_par[4], est_par[5], est_par[6]), nsim = 50)
             set.seed(3)
             g3_sim_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(est_par[3], est_par[4], est_par[5], est_par[6]), nsim = 50)
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- stpm(data.frame(x = c(g1_sim_ppp[[si]]$df$x, g2_sim_ppp[[si]]$df$x, g3_sim_ppp[[si]]$df$x),
                                                y = c(g1_sim_ppp[[si]]$df$y, g2_sim_ppp[[si]]$df$y, g3_sim_ppp[[si]]$df$y),
                                                t = c(g1_sim_ppp[[si]]$df$t, g2_sim_ppp[[si]]$df$t, g3_sim_ppp[[si]]$df$t),
                                                m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]]$df)), rep("b", nrow(g2_sim_ppp[[si]]$df)), rep("c", nrow(g3_sim_ppp[[si]]$df))))))
             }
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$df$m1
               linear_predictor_sim <- param_vec[marks_vec_sim] + est_par[4] * sim_ppp[[j]]$df$x + est_par[5] * sim_ppp[[j]]$df$y + est_par[6] * sim_ppp[[j]]$df$t
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           },
           
           "8" = {
             est <- stppm_prova_local(X = ppp_true, formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = comb$mult, ncube = comb$ncube, grid = comb$grid, marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)
             est_par <- c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]), 
                          est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
                          est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[8]),
                          est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), 
                          est_param6 = as.numeric(est$IntCoefs[4]), est_param7 = as.numeric(est$IntCoefs[5]))
             Khat_data <- Khat_spatial.3D(stp(ppp_true$df[,c(1,2,3)]), lambda = est$l, correction = "translate")
             set.seed(1)
             g1_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[1], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             set.seed(2)
             g2_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[2], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             set.seed(3)
             g3_sim_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                    par = c(est_par[3], est_par[4], est_par[5], est_par[6], est_par[7]), nsim = 50, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
             sim_ppp <- setNames(replicate(50, data.frame()), paste0("X", 1:50))
             for(si in 1:50) {
               sim_ppp[[si]] <- data.frame(x = c(g1_sim_ppp[[si]]$x, g2_sim_ppp[[si]]$x, g3_sim_ppp[[si]]$x),
                                           y = c(g1_sim_ppp[[si]]$y, g2_sim_ppp[[si]]$y, g3_sim_ppp[[si]]$y),
                                           t = c(g1_sim_ppp[[si]]$t, g2_sim_ppp[[si]]$t, g3_sim_ppp[[si]]$t),
                                           cov1 = c(g1_sim_ppp[[si]]$V1, g2_sim_ppp[[si]]$V1, g3_sim_ppp[[si]]$V1),
                                           m1 = as.factor(c(rep("a", nrow(g1_sim_ppp[[si]])), rep("b", nrow(g2_sim_ppp[[si]])), rep("c", nrow(g3_sim_ppp[[si]])))))
             }
             sim_ppp <- Filter(function(df) all(complete.cases(df)), sim_ppp)
             param_vec <- c(a = est_par[1], b = est_par[2], c = est_par[3])
             Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_data$dist), nrow = length(sim_ppp)))
             for(j in 1:length(sim_ppp)){
               marks_vec_sim <- sim_ppp[[j]]$m1
               linear_predictor_sim <- param_vec[marks_vec_sim] + est_par[4] * sim_ppp[[j]]$x + est_par[5] * sim_ppp[[j]]$y + est_par[6] * sim_ppp[[j]]$t + est_par[7] * sim_ppp[[j]]$cov1
               pred_int <- exp(linear_predictor_sim)
               Khat_sim_j <- Khat_spatial.3D(stp(sim_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data$dist)
               Khat_sim[j, ] <- Khat_sim_j$Khat
             }
             curve_set <- create_curve_set(list(r = Khat_data$dist, obs = Khat_data$Khat, sim_m = t(Khat_sim)))
             global_env <- global_envelope_test(curve_set, alternative = "two.sided")
             pvalue <- attr(global_env, "p")
             pvalue
           }
    )
  }, error = function(e) NA_real_)
  
  p() 
  return(result)
}

with_progress({
  p <- progressor(steps = n_comb_total)
  results_n500_pvalue <- future_lapply(1:n_comb_total, calc_pvalue_n500, p = p)
})
dataset_n500_pvalue$pvalue <- unlist(results_n500_pvalue)
save(dataset_n500_pvalue, file = "dataset_n500_pvalue_parallel.RData")



#### ANALISI RISULTATI VALIDAZIONE LINEE GUIDA ####

dataset_n100_pvalue_sig <- dataset_n100_pvalue[which(dataset_n100_pvalue$pvalue < 0.05),]
dataset_n100_pvalue_NA <- dataset_n100_pvalue[which(is.na(dataset_n100_pvalue$pvalue)),]
table(dataset_n100_pvalue_NA$grid)

dataset_n250_pvalue_NA <- dataset_n250_pvalue[which(is.na(dataset_n250_pvalue$pvalue)),]
table(dataset_n250_pvalue_NA$grid)

dataset_n500_pvalue_NA <- dataset_n100_pvalue[which(is.na(dataset_n500_pvalue$pvalue)),]
table(dataset_n500_pvalue_NA$grid)

# --- 1️⃣ Unisci tutti i dataset ---
dataset_all_validation <- bind_rows(
  dataset_n100_pvalue %>% mutate(n = 100),
  dataset_n250_pvalue %>% mutate(n = 250),
  dataset_n500_pvalue %>% mutate(n = 500)
)

# --- 2️⃣ Crea categorie di p-value e calcola percentuali per ciascun set ---
dataset_summary_all_validation <- dataset_all_validation %>%
  mutate(
    p_cat = case_when(
      pvalue < 0.01 ~ "<0.01",
      pvalue < 0.05 ~ "0.01–0.05",
      pvalue < 0.1  ~ "0.05–0.1",
      TRUE ~ ">0.1"
    ),
    p_cat = factor(p_cat, levels = c("<0.01", "0.01–0.05", "0.05–0.1", ">0.1"))
  ) %>%
  group_by(grid, n, p_cat) %>%
  summarise(freq = n(), .groups = "drop_last") %>%
  mutate(prop = 100 * freq / sum(freq))

# --- 3️⃣ Grafico: istogramma con 4 categorie di p-value ---
ggplot(dataset_summary_all_validation, aes(x = p_cat, y = prop, fill = p_cat)) +
  geom_col(color = "black", width = 0.7) +
  facet_grid(grid ~ n, labeller = label_both) +
  scale_fill_manual(
    values = c("<0.01" = "red",
               "0.01–0.05" = "orange",
               "0.05–0.1" = "gold",
               ">0.1" = "white")
  ) +
  labs(
    x = "p-value",
    y = "Percentage",
    fill = "p-value class",
    title = ""
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )


table_summary <- dataset_all_validation %>%
  group_by(grid, n) %>%
  summarise(
    n_comb = n_distinct(paste(mult, ncube, sep = "-")), # numero combinazioni mult-ncube
    perc_p_gt_0.05 = 100 * mean(pvalue > 0.05, na.rm = T),
    perc_p_gt_0.1 = 100 * mean(pvalue > 0.1, na.rm = T),
    .groups = "drop"
  ) %>%
  arrange(n, grid)

