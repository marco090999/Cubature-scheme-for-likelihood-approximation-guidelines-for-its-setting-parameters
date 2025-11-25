

######################################
#### CODES SIMULATION STUDY - MSE ####
######################################

library(gridExtra)
library(furrr)
library(tibble)
library(dplyr)
library(purrr)
library(progressr)
library(parallel)
library(pbapply)
library(spatstat)
library(stpp) 
library(stopp)
library(mgcv)

n_cores <- 64
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(spatstat)
  library(stpp) 
  library(stopp)
  library(mgcv)
})

#################################################
#### SCENARIO 1: POISSON INHOM CON X N = 100 ####
#################################################

set.seed(2)
inh_scen1_n100 <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(3.5, 2), nsim = 100)

n <- NULL
for(i in 1:length(inh_scen1_n100)) {
  n[i] <- nrow(inh_scen1_n100[[i]]$df)
}
mean(n)

true_params_scen1_n100 <- c(3.5, 2)
mult_param_scen1_n100 <- seq(1, 50, by = 1)
ncube_param_scen1_n100 <- seq(1, 50, by = 1)
grid_param_scen1_n100 <- c(FALSE, TRUE)
df_param_scen1_n100 <- expand.grid(mult = mult_param_scen1_n100, ncube = ncube_param_scen1_n100, grid = grid_param_scen1_n100)

n_settings_scen1_n100 <- nrow(df_param_scen1_n100)
n_params_scen1_n100 <- length(true_params_scen1_n100)

results_list_scen1_n100 <- vector("list", length = length(inh_scen1_n100))

for (i in seq_along(inh_scen1_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen1_n100[[i]]
  
  export_list_scen1_n100 <- c("process_chunk","sim_process", "df_param_scen1_n100", "inh_scen1_n100", "true_params_scen1_n100", functions_in_workspace)
  clusterExport(cl, varlist = export_list_scen1_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen1_n100, function(j) {
    mult_val <- df_param_scen1_n100$mult[j]
    ncube_val <- df_param_scen1_n100$ncube[j]
    grid_val <- df_param_scen1_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x, mult = mult_val, ncube = ncube_val, grid = grid_val, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2])))
  })
  
  res_df <- cbind(df_param_scen1_n100, t(res_matrix))
  results_list_scen1_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen1_n100, file = "results_list_scen1_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 1 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen1_n100), elapsed_time))
}


all_results_scen1_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen1_n100, seq_along(results_list_scen1_n100)))

summary_stats_scen1_n100 <- all_results_scen1_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    bias_param1 = mean_param1 - true_params_scen1_n100[1],
    bias_param2 = mean_param2 - true_params_scen1_n100[2],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_combined = (mse_param1 + mse_param2) / 2,
    .groups = 'drop'
  )



#################################################
#### SCENARIO 1: POISSON INHOM CON X N = 250 ####
#################################################

set.seed(2)
inh_scen1_n250 <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(4.3, 2), nsim = 100)

n <- NULL
for(i in 1:length(inh_scen1_n250)) {
  n[i] <- nrow(inh_scen1_n250[[i]]$df)
}
mean(n)

true_params_scen1_n250 <- c(4.3, 2)
mult_param_scen1_n250 <- seq(1, 50, by = 1)
ncube_param_scen1_n250 <- seq(1, 50, by = 1)
grid_param_scen1_n250 <- c(FALSE, TRUE)
df_param_scen1_n250 <- expand.grid(mult = mult_param_scen1_n250, ncube = ncube_param_scen1_n250, grid = grid_param_scen1_n250)

n_settings_scen1_n250 <- nrow(df_param_scen1_n250)
n_params_scen1_n250 <- length(true_params_scen1_n250)

results_list_scen1_n250 <- vector("list", length = length(inh_scen1_n250))

for (i in seq_along(inh_scen1_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen1_n250[[i]]
  
  export_list_scen1_n250 <- c("process_chunk","sim_process", "df_param_scen1_n250", "inh_scen1_n250", "true_params_scen1_n250", functions_in_workspace)
  clusterExport(cl, varlist = export_list_scen1_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen1_n250, function(j) {
    mult_val <- df_param_scen1_n250$mult[j]
    ncube_val <- df_param_scen1_n250$ncube[j]
    grid_val <- df_param_scen1_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x, mult = mult_val, ncube = ncube_val, grid = grid_val, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2])))
  })
  
  res_df <- cbind(df_param_scen1_n250, t(res_matrix))
  results_list_scen1_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen1_n250, file = "results_list_scen1_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen1_n250), elapsed_time))
}


all_results_scen1_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen1_n250, seq_along(results_list_scen1_n250)))

summary_stats_scen1_n250 <- all_results_scen1_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    bias_param1 = mean_param1 - true_params_scen1_n250[1],
    bias_param2 = mean_param2 - true_params_scen1_n250[2],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_combined = (mse_param1 + mse_param2) / 2,
    .groups = 'drop'
  )



#################################################
#### SCENARIO 1: POISSON INHOM CON X N = 500 ####
#################################################

set.seed(2)
inh_scen1_n500 <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(5, 2), nsim = 100)

n <- NULL
for(i in 1:length(inh_scen1_n500)) {
  n[i] <- nrow(inh_scen1_n500[[i]]$df)
}
mean(n)

true_params_scen1_n500 <- c(5, 2)
mult_param_scen1_n500 <- seq(1, 50, by = 1)
ncube_param_scen1_n500 <- seq(1, 50, by = 1)
grid_param_scen1_n500 <- c(FALSE, TRUE)
df_param_scen1_n500 <- expand.grid(mult = mult_param_scen1_n500, ncube = ncube_param_scen1_n500, grid = grid_param_scen1_n500)

n_settings_scen1_n500 <- nrow(df_param_scen1_n500)
n_params_scen1_n500 <- length(true_params_scen1_n500)

results_list_scen1_n500 <- vector("list", length = length(inh_scen1_n500))

for (i in seq_along(inh_scen1_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen1_n500[[i]]
  
  export_list_scen1_n500 <- c("process_chunk","sim_process", "df_param_scen1_n500", "inh_scen1_n500", "true_params_scen1_n500", functions_in_workspace)
  clusterExport(cl, varlist = export_list_scen1_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen1_n500, function(j) {
    mult_val <- df_param_scen1_n500$mult[j]
    ncube_val <- df_param_scen1_n500$ncube[j]
    grid_val <- df_param_scen1_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x, mult = mult_val, ncube = ncube_val, grid = grid_val, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2])))
  })
  
  res_df <- cbind(df_param_scen1_n500, t(res_matrix))
  results_list_scen1_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen1_n500, file = "results_list_scen1_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen1_n500), elapsed_time))
}


all_results_scen1_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen1_n500, seq_along(results_list_scen1_n500)))

summary_stats_scen1_n500 <- all_results_scen1_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    bias_param1 = mean_param1 - true_params_scen1_n500[1],
    bias_param2 = mean_param2 - true_params_scen1_n500[2],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_combined = (mse_param1 + mse_param2) / 2,
    .groups = 'drop'
  )



###################################################
#### SCENARIO 2: POISSON INHOM CON XYT N = 500 ####
###################################################

set.seed(2)
inh_scen2_n500 <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(2.5, 2, 2, 2), nsim = 100)

n <- NULL
for(i in 1:length(inh_scen2_n500)) {
  n[i] <- nrow(inh_scen2_n500[[i]]$df)
}
mean(n)

true_params_scen2_n500 <- c(2.5, 2, 2, 2)
mult_param_scen2_n500 <- seq(1, 50, by = 1)
ncube_param_scen2_n500 <- seq(1, 50, by = 1)
grid_param_scen2_n500 <- c(FALSE, TRUE)
df_param_scen2_n500 <- expand.grid(mult = mult_param_scen2_n500, ncube = ncube_param_scen2_n500, grid = grid_param_scen2_n500)

n_settings_scen2_n500 <- nrow(df_param_scen2_n500)
n_params_scen2_n500 <- length(true_params_scen2_n500)

results_list_scen2_n500 <- vector("list", length = length(inh_scen2_n500))

for (i in seq_along(inh_scen2_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen2_n500[[i]]
  
  export_list_scen2_n500 <- c("process_chunk","sim_process", "df_param_scen2_n500", "inh_scen2_n500", "true_params_scen2_n500", functions_in_workspace)
  clusterExport(cl, varlist = export_list_scen2_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen2_n500, function(j) {
    mult_val <- df_param_scen2_n500$mult[j]
    ncube_val <- df_param_scen2_n500$ncube[j]
    grid_val <- df_param_scen2_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t, mult = mult_val, ncube = ncube_val, grid = grid_val, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen2_n500, t(res_matrix))
  results_list_scen2_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen2_n500, file = "results_list_scen2_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen2_n500), elapsed_time))
}


all_results_scen2_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen2_n500, seq_along(results_list_scen2_n500)))

summary_stats_scen2_n500 <- all_results_scen2_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    bias_param1 = mean_param1 - true_params_scen2_n500[1],
    bias_param2 = mean_param2 - true_params_scen2_n500[2],
    bias_param3 = mean_param3 - true_params_scen2_n500[3],
    bias_param4 = mean_param4 - true_params_scen2_n500[4],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4) / 4,
    .groups = 'drop'
  )



###################################################
#### SCENARIO 2: POISSON INHOM CON XYT N = 250 ####
###################################################

set.seed(2)
inh_scen2_n250 <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(1.8, 2, 2, 2), nsim = 100)

n <- NULL
for(i in 1:length(inh_scen2_n250)) {
  n[i] <- nrow(inh_scen2_n250[[i]]$df)
}
mean(n)

true_params_scen2_n250 <- c(1.8, 2, 2, 2)
mult_param_scen2_n250 <- seq(1, 50, by = 1)
ncube_param_scen2_n250 <- seq(1, 50, by = 1)
grid_param_scen2_n250 <- c(FALSE, TRUE)
df_param_scen2_n250 <- expand.grid(mult = mult_param_scen2_n250, ncube = ncube_param_scen2_n250, grid = grid_param_scen2_n250)

n_settings_scen2_n250 <- nrow(df_param_scen2_n250)
n_params_scen2_n250 <- length(true_params_scen2_n250)

results_list_scen2_n250 <- vector("list", length = length(inh_scen2_n250))

for (i in seq_along(inh_scen2_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen2_n250[[i]]
  
  export_list_scen2_n250 <- c("process_chunk","sim_process", "df_param_scen2_n250", "inh_scen2_n250", "true_params_scen2_n250", functions_in_workspace)
  clusterExport(cl, varlist = export_list_scen2_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen2_n250, function(j) {
    mult_val <- df_param_scen2_n250$mult[j]
    ncube_val <- df_param_scen2_n250$ncube[j]
    grid_val <- df_param_scen2_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t, mult = mult_val, ncube = ncube_val, grid = grid_val, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen2_n250, t(res_matrix))
  results_list_scen2_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen2_n250, file = "results_list_scen2_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen2_n250), elapsed_time))
}


all_results_scen2_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen2_n250, seq_along(results_list_scen2_n250)))

summary_stats_scen2_n250 <- all_results_scen2_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    bias_param1 = mean_param1 - true_params_scen2_n250[1],
    bias_param2 = mean_param2 - true_params_scen2_n250[2],
    bias_param3 = mean_param3 - true_params_scen2_n250[3],
    bias_param4 = mean_param4 - true_params_scen2_n250[4],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4) / 4,
    .groups = 'drop'
  )



###################################################
#### SCENARIO 2: POISSON INHOM CON XYT N = 100 ####
###################################################

set.seed(2)
inh_scen2_n100 <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(1, 2, 2, 2), nsim = 100)

n <- NULL
for(i in 1:length(inh_scen2_n100)) {
  n[i] <- nrow(inh_scen2_n100[[i]]$df)
}
mean(n)

true_params_scen2_n100 <- c(1, 2, 2, 2)
mult_param_scen2_n100 <- seq(1, 50, by = 1)
ncube_param_scen2_n100 <- seq(1, 50, by = 1)
grid_param_scen2_n100 <- c(FALSE, TRUE)
df_param_scen2_n100 <- expand.grid(mult = mult_param_scen2_n100, ncube = ncube_param_scen2_n100, grid = grid_param_scen2_n100)

n_settings_scen2_n100 <- nrow(df_param_scen2_n100)
n_params_scen2_n100 <- length(true_params_scen2_n100)

results_list_scen2_n100 <- vector("list", length = length(inh_scen2_n100))

for (i in seq_along(inh_scen2_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen2_n100[[i]]
  
  export_list_scen2_n100 <- c("process_chunk","sim_process", "df_param_scen2_n100", "inh_scen2_n100", "true_params_scen2_n100", functions_in_workspace)
  clusterExport(cl, varlist = export_list_scen2_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen2_n100, function(j) {
    mult_val <- df_param_scen2_n100$mult[j]
    ncube_val <- df_param_scen2_n100$ncube[j]
    grid_val <- df_param_scen2_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t, mult = mult_val, ncube = ncube_val, grid = grid_val, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen2_n100, t(res_matrix))
  results_list_scen2_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen2_n100, file = "results_list_scen2_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen2_n100), elapsed_time))
}


all_results_scen2_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen2_n100, seq_along(results_list_scen2_n100)))

summary_stats_scen2_n100 <- all_results_scen2_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    bias_param1 = mean_param1 - true_params_scen2_n100[1],
    bias_param2 = mean_param2 - true_params_scen2_n100[2],
    bias_param3 = mean_param3 - true_params_scen2_n100[3],
    bias_param4 = mean_param4 - true_params_scen2_n100[4],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4) / 4,
    .groups = 'drop'
  )



##################################################
#### SCENARIO 3: POISSON INHOM CON Z1 N = 100 ####
##################################################

sim.cov <- function(n.obs = 200, minX = 0, maxX = 1, minY = 0, maxY = 1, minT = 0, maxT = 1) {
  x <- runif(n.obs, min =  minX, max = maxX)
  y <- runif(n.obs, min =  minY, max = maxY)
  t <- runif(n.obs, min =  minT, max = maxT)
  
  mean1 <- 2*x + 2*y + 2*t
  mean2 <- 2*x + 2*y - 2*t
  mean3 <- 2*x - 2*y - 2*t
  
  cov1 <- rnorm(n.obs, mean = mean1, sd = 1)
  cov2 <- rnorm(n.obs, mean = mean2, sd = 1)
  cov3 <- rnorm(n.obs, mean = mean3, sd = 1)
  
  cov.list <- list(
    cov1 = data.frame(x = x, y = y, t = t, cov1 = cov1),
    cov2 = data.frame(x = x, y = y, t = t, cov2 = cov2),
    cov3 = data.frame(x = x, y = y, t = t, cov3 = cov3)
  )
  
  return(cov.list)
}

set.seed(2)
df_covs = sim.cov()

df_covs[[1]]$cov1 <- scale_to_range(df_covs[[1]]$cov1, 0, 1)
df_covs[[2]]$cov2 <- scale_to_range(df_covs[[2]]$cov2, 0, 1)
df_covs[[3]]$cov3 <- scale_to_range(df_covs[[3]]$cov3, 0, 1)

plot3D::scatter3D(df_covs[[1]]$x, df_covs[[1]]$y, df_covs[[1]]$t, colvar = df_covs[[1]]$cov)
plot3D::scatter3D(df_covs[[2]]$x, df_covs[[2]]$y, df_covs[[2]]$t, colvar = df_covs[[2]]$cov)
plot3D::scatter3D(df_covs[[3]]$x, df_covs[[3]]$y, df_covs[[3]]$t, colvar = df_covs[[3]]$cov)

inh_scen3_n100 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])}, par = c(3.6, 2), 
                                        covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen3_n100)) {
  n[i] <- nrow(inh_scen3_n100[[i]]$df)
}
mean(n)

true_params_scen3_n100 <- c(3.6, 2)
mult_param_scen3_n100 <- seq(1, 50, by = 1)
ncube_param_scen3_n100 <- seq(1, 50, by = 1)
grid_param_scen3_n100 <- c(FALSE, TRUE)
df_param_scen3_n100 <- expand.grid(mult = mult_param_scen3_n100, ncube = ncube_param_scen3_n100, grid = grid_param_scen3_n100)

n_settings_scen3_n100 <- nrow(df_param_scen3_n100)
n_params_scen3_n100 <- length(true_params_scen3_n100)

results_list_scen3_n100 <- vector("list", length = length(inh_scen3_n100))

for (i in seq_along(inh_scen3_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen3_n100[[i]]
  
  export_list_scen3_n100 <- c("sim_process", "df_param_scen3_n100", "inh_scen3_n100", "true_params_scen3_n100", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen3_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen3_n100, function(j) {
    mult_val <- df_param_scen3_n100$mult[j]
    ncube_val <- df_param_scen3_n100$ncube[j]
    grid_val <- df_param_scen3_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ cov1, mult = mult_val, ncube = ncube_val, grid = grid_val, spatial.cov = TRUE,
                             covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2])))
  })
  
  res_df <- cbind(df_param_scen3_n100, t(res_matrix))
  results_list_scen3_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen3_n100, file = "results_list_scen3_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 3 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen3_n100), elapsed_time))
}


all_results_scen3_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen3_n100, seq_along(results_list_scen3_n100)))

summary_stats_scen3_n100 <- all_results_scen3_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    bias_param1 = mean_param1 - true_params_scen3_n100[1],
    bias_param2 = mean_param2 - true_params_scen3_n100[2],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_combined = (mse_param1 + mse_param2) / 2,
    .groups = 'drop'
  )



##################################################
#### SCENARIO 3: POISSON INHOM CON Z1 N = 250 ####
##################################################

inh_scen3_n250 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])}, par = c(4.5, 2), 
                                        covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen3_n250)) {
  n[i] <- nrow(inh_scen3_n250[[i]]$df)
}
mean(n)

true_params_scen3_n250 <- c(4.5, 2)
mult_param_scen3_n250 <- seq(1, 50, by = 1)
ncube_param_scen3_n250 <- seq(1, 50, by = 1)
grid_param_scen3_n250 <- c(FALSE, TRUE)
df_param_scen3_n250 <- expand.grid(mult = mult_param_scen3_n250, ncube = ncube_param_scen3_n250, grid = grid_param_scen3_n250)

n_settings_scen3_n250 <- nrow(df_param_scen3_n250)
n_params_scen3_n250 <- length(true_params_scen3_n250)

results_list_scen3_n250 <- vector("list", length = length(inh_scen3_n250))

for (i in seq_along(inh_scen3_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen3_n250[[i]]
  
  export_list_scen3_n250 <- c("sim_process", "df_param_scen3_n250", "inh_scen3_n250", "true_params_scen3_n250", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen3_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen3_n250, function(j) {
    mult_val <- df_param_scen3_n250$mult[j]
    ncube_val <- df_param_scen3_n250$ncube[j]
    grid_val <- df_param_scen3_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ cov1, mult = mult_val, ncube = ncube_val, grid = grid_val, spatial.cov = TRUE,
                             covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2])))
  })
  
  res_df <- cbind(df_param_scen3_n250, t(res_matrix))
  results_list_scen3_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen3_n250, file = "results_list_scen3_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 3 n = 250: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen3_n250), elapsed_time))
}


all_results_scen3_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen3_n250, seq_along(results_list_scen3_n250)))

summary_stats_scen3_n250 <- all_results_scen3_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    bias_param1 = mean_param1 - true_params_scen3_n250[1],
    bias_param2 = mean_param2 - true_params_scen3_n250[2],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_combined = (mse_param1 + mse_param2) / 2,
    .groups = 'drop'
  )



##################################################
#### SCENARIO 3: POISSON INHOM CON Z1 N = 500 ####
##################################################

inh_scen3_n500 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])}, par = c(5.3, 2), 
                                        covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen3_n500)) {
  n[i] <- nrow(inh_scen3_n500[[i]]$df)
}
mean(n)

true_params_scen3_n500 <- c(5.3, 2)
mult_param_scen3_n500 <- seq(1, 50, by = 1)
ncube_param_scen3_n500 <- seq(1, 50, by = 1)
grid_param_scen3_n500 <- c(FALSE, TRUE)
df_param_scen3_n500 <- expand.grid(mult = mult_param_scen3_n500, ncube = ncube_param_scen3_n500, grid = grid_param_scen3_n500)

n_settings_scen3_n500 <- nrow(df_param_scen3_n500)
n_params_scen3_n500 <- length(true_params_scen3_n500)

results_list_scen3_n500 <- vector("list", length = length(inh_scen3_n500))

for (i in 55:100) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen3_n500[[i]]
  
  export_list_scen3_n500 <- c("sim_process", "df_param_scen3_n500", "inh_scen3_n500", "true_params_scen3_n500", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen3_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen3_n500, function(j) {
    mult_val <- df_param_scen3_n500$mult[j]
    ncube_val <- df_param_scen3_n500$ncube[j]
    grid_val <- df_param_scen3_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ cov1, mult = mult_val, ncube = ncube_val, grid = grid_val, spatial.cov = TRUE,
                             covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2])))
  })
  
  res_df <- cbind(df_param_scen3_n500, t(res_matrix))
  results_list_scen3_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen3_n500, file = "results_list_scen3_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 3 n = 500: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen3_n500), elapsed_time))
}


all_results_scen3_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen3_n500, seq_along(results_list_scen3_n500)))

summary_stats_scen3_n500 <- all_results_scen3_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    bias_param1 = mean_param1 - true_params_scen3_n500[1],
    bias_param2 = mean_param2 - true_params_scen3_n500[2],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_combined = (mse_param1 + mse_param2) / 2,
    .groups = 'drop'
  )



########################################################
#### SCENARIO 4: POISSON INHOM CON Z1-Z2-Z3 N = 100 ####
########################################################

inh_scen4_n100 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, 
                                        par = c(1.4, 2, 2, 2), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen4_n100)) {
  n[i] <- nrow(inh_scen4_n100[[i]]$df)
}
mean(n)

true_params_scen4_n100 <- c(1.4, 2, 2, 2)
mult_param_scen4_n100 <- seq(1, 50, by = 1)
ncube_param_scen4_n100 <- seq(1, 50, by = 1)
grid_param_scen4_n100 <- c(FALSE, TRUE)
df_param_scen4_n100 <- expand.grid(mult = mult_param_scen4_n100, ncube = ncube_param_scen4_n100, grid = grid_param_scen4_n100)

n_settings_scen4_n100 <- nrow(df_param_scen4_n100)
n_params_scen4_n100 <- length(true_params_scen4_n100)

results_list_scen4_n100 <- vector("list", length = length(inh_scen4_n100))

for (i in seq_along(inh_scen4_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen4_n100[[i]]
  
  export_list_scen4_n100 <- c("sim_process", "df_param_scen4_n100", "inh_scen4_n100", "true_params_scen4_n100", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen4_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen4_n100, function(j) {
    mult_val <- df_param_scen4_n100$mult[j]
    ncube_val <- df_param_scen4_n100$ncube[j]
    grid_val <- df_param_scen4_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ cov1 + cov2 + cov3, mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen4_n100, t(res_matrix))
  results_list_scen4_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen4_n100, file = "results_list_scen4_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 4 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen4_n100), elapsed_time))
}


all_results_scen4_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen4_n100, seq_along(results_list_scen4_n100)))

summary_stats_scen4_n100 <- all_results_scen4_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    bias_param1 = mean_param1 - true_params_scen4_n100[1],
    bias_param2 = mean_param2 - true_params_scen4_n100[2],
    bias_param3 = mean_param3 - true_params_scen4_n100[3],
    bias_param4 = mean_param4 - true_params_scen4_n100[4],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4) / 4,
    .groups = 'drop'
  )



########################################################
#### SCENARIO 4: POISSON INHOM CON Z1-Z2-Z3 N = 250 ####
########################################################

inh_scen4_n250 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, 
                                        par = c(2.2, 2, 2, 2), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen4_n250)) {
  n[i] <- nrow(inh_scen4_n250[[i]]$df)
}
mean(n)

true_params_scen4_n250 <- c(2.2, 2, 2, 2)
mult_param_scen4_n250 <- seq(1, 50, by = 1)
ncube_param_scen4_n250 <- seq(1, 50, by = 1)
grid_param_scen4_n250 <- c(FALSE, TRUE)
df_param_scen4_n250 <- expand.grid(mult = mult_param_scen4_n250, ncube = ncube_param_scen4_n250, grid = grid_param_scen4_n250)

n_settings_scen4_n250 <- nrow(df_param_scen4_n250)
n_params_scen4_n250 <- length(true_params_scen4_n250)

results_list_scen4_n250 <- vector("list", length = length(inh_scen4_n250))

for (i in 84:100) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen4_n250[[i]]
  
  export_list_scen4_n250 <- c("sim_process", "df_param_scen4_n250", "inh_scen4_n250", "true_params_scen4_n250", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen4_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen4_n250, function(j) {
    mult_val <- df_param_scen4_n250$mult[j]
    ncube_val <- df_param_scen4_n250$ncube[j]
    grid_val <- df_param_scen4_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ cov1 + cov2 + cov3, mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen4_n250, t(res_matrix))
  results_list_scen4_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen4_n250, file = "results_list_scen4_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 4 n = 250: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen4_n250), elapsed_time))
}


all_results_scen4_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen4_n250, seq_along(results_list_scen4_n250)))

summary_stats_scen4_n250 <- all_results_scen4_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    bias_param1 = mean_param1 - true_params_scen4_n250[1],
    bias_param2 = mean_param2 - true_params_scen4_n250[2],
    bias_param3 = mean_param3 - true_params_scen4_n250[3],
    bias_param4 = mean_param4 - true_params_scen4_n250[4],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4) / 4,
    .groups = 'drop'
  )



########################################################
#### SCENARIO 4: POISSON INHOM CON Z1-Z2-Z3 N = 500 ####
########################################################

inh_scen4_n500 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, 
                                        par = c(2.9, 2, 2, 2), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen4_n500)) {
  n[i] <- nrow(inh_scen4_n500[[i]]$df)
}
mean(n)

true_params_scen4_n500 <- c(2.9, 2, 2, 2)
mult_param_scen4_n500 <- seq(1, 50, by = 1)
ncube_param_scen4_n500 <- seq(1, 50, by = 1)
grid_param_scen4_n500 <- c(FALSE, TRUE)
df_param_scen4_n500 <- expand.grid(mult = mult_param_scen4_n500, ncube = ncube_param_scen4_n500, grid = grid_param_scen4_n500)

n_settings_scen4_n500 <- nrow(df_param_scen4_n500)
n_params_scen4_n500 <- length(true_params_scen4_n500)

results_list_scen4_n500 <- vector("list", length = length(inh_scen4_n500))

for (i in seq_along(inh_scen4_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen4_n500[[i]]
  
  export_list_scen4_n500 <- c("sim_process", "df_param_scen4_n500", "inh_scen4_n500", "true_params_scen4_n500", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen4_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen4_n500, function(j) {
    mult_val <- df_param_scen4_n500$mult[j]
    ncube_val <- df_param_scen4_n500$ncube[j]
    grid_val <- df_param_scen4_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ cov1 + cov2 + cov3, mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen4_n500, t(res_matrix))
  results_list_scen4_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen4_n500, file = "results_list_scen4_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 4 n = 500: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen4_n500), elapsed_time))
}


all_results_scen4_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen4_n500, seq_along(results_list_scen4_n500)))

summary_stats_scen4_n500 <- all_results_scen4_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    bias_param1 = mean_param1 - true_params_scen4_n500[1],
    bias_param2 = mean_param2 - true_params_scen4_n500[2],
    bias_param3 = mean_param3 - true_params_scen4_n500[3],
    bias_param4 = mean_param4 - true_params_scen4_n500[4],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4) / 4,
    .groups = 'drop'
  )



############################################################
#### SCENARIO 5: POISSON INHOM CON xyt Z1-Z2-Z3 N = 100 ####
############################################################

inh_scen5_n100 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])}, 
                                        par = c(3, 2, 2, 2, -1, -1, -1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen5_n100)) {
  n[i] <- nrow(inh_scen5_n100[[i]]$df)
}
mean(n)

true_params_scen5_n100 <- c(3, 2, 2, 2, -1, -1, -1)
mult_param_scen5_n100 <- seq(1, 50, by = 1)
ncube_param_scen5_n100 <- seq(1, 50, by = 1)
grid_param_scen5_n100 <- c(FALSE, TRUE)
df_param_scen5_n100 <- expand.grid(mult = mult_param_scen5_n100, ncube = ncube_param_scen5_n100, grid = grid_param_scen5_n100)

n_settings_scen5_n100 <- nrow(df_param_scen5_n100)
n_params_scen5_n100 <- length(true_params_scen5_n100)

results_list_scen5_n100 <- vector("list", length = length(inh_scen5_n100))

for (i in seq_along(inh_scen5_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen5_n100[[i]]
  
  export_list_scen5_n100 <- c("sim_process", "df_param_scen5_n100", "inh_scen5_n100", "true_params_scen5_n100", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen5_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen5_n100, function(j) {
    mult_val <- df_param_scen5_n100$mult[j]
    ncube_val <- df_param_scen5_n100$ncube[j]
    grid_val <- df_param_scen5_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + cov1 + cov2 + cov3, mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]),
             est_param5 = as.numeric(est$IntCoefs[5]), est_param6 = as.numeric(est$IntCoefs[6]),
             est_param7 = as.numeric(est$IntCoefs[7])))
  })
  
  res_df <- cbind(df_param_scen5_n100, t(res_matrix))
  results_list_scen5_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen5_n100, file = "results_list_scen5_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 5 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen5_n100), elapsed_time))
}


all_results_scen5_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen5_n100, seq_along(results_list_scen5_n100)))

summary_stats_scen5_n100 <- all_results_scen5_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    mean_param7 = mean(est_param7),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    sd_param7 = sd(est_param7),
    bias_param1 = mean_param1 - true_params_scen5_n100[1],
    bias_param2 = mean_param2 - true_params_scen5_n100[2],
    bias_param3 = mean_param3 - true_params_scen5_n100[3],
    bias_param4 = mean_param4 - true_params_scen5_n100[4],
    bias_param5 = mean_param5 - true_params_scen5_n100[5],
    bias_param6 = mean_param6 - true_params_scen5_n100[6],
    bias_param7 = mean_param7 - true_params_scen5_n100[7],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_param7 = sd_param7^2 + bias_param7^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6 + mse_param7) / 7,
    .groups = 'drop'
  )



############################################################
#### SCENARIO 5: POISSON INHOM CON xyt Z1-Z2-Z3 N = 250 ####
############################################################

inh_scen5_n250 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])}, 
                                        par = c(3.7, 2, 2, 2, -1, -1, -1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen5_n250)) {
  n[i] <- nrow(inh_scen5_n250[[i]]$df)
}
mean(n)

true_params_scen5_n250 <- c(3.7, 2, 2, 2, -1, -1, -1)
mult_param_scen5_n250 <- seq(1, 50, by = 1)
ncube_param_scen5_n250 <- seq(1, 50, by = 1)
grid_param_scen5_n250 <- c(FALSE, TRUE)
df_param_scen5_n250 <- expand.grid(mult = mult_param_scen5_n250, ncube = ncube_param_scen5_n250, grid = grid_param_scen5_n250)

n_settings_scen5_n250 <- nrow(df_param_scen5_n250)
n_params_scen5_n250 <- length(true_params_scen5_n250)

results_list_scen5_n250 <- vector("list", length = length(inh_scen5_n250))

for (i in seq_along(inh_scen5_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen5_n250[[i]]
  
  export_list_scen5_n250 <- c("sim_process", "df_param_scen5_n250", "inh_scen5_n250", "true_params_scen5_n250", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen5_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen5_n250, function(j) {
    mult_val <- df_param_scen5_n250$mult[j]
    ncube_val <- df_param_scen5_n250$ncube[j]
    grid_val <- df_param_scen5_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + cov1 + cov2 + cov3, mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]),
             est_param5 = as.numeric(est$IntCoefs[5]), est_param6 = as.numeric(est$IntCoefs[6]),
             est_param7 = as.numeric(est$IntCoefs[7])))
  })
  
  res_df <- cbind(df_param_scen5_n250, t(res_matrix))
  results_list_scen5_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen5_n250, file = "results_list_scen5_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 5 n = 250: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen5_n250), elapsed_time))
}


all_results_scen5_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen5_n250, seq_along(results_list_scen5_n250)))

summary_stats_scen5_n250 <- all_results_scen5_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    mean_param7 = mean(est_param7),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    sd_param7 = sd(est_param7),
    bias_param1 = mean_param1 - true_params_scen5_n250[1],
    bias_param2 = mean_param2 - true_params_scen5_n250[2],
    bias_param3 = mean_param3 - true_params_scen5_n250[3],
    bias_param4 = mean_param4 - true_params_scen5_n250[4],
    bias_param5 = mean_param5 - true_params_scen5_n250[5],
    bias_param6 = mean_param6 - true_params_scen5_n250[6],
    bias_param7 = mean_param7 - true_params_scen5_n250[7],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_param7 = sd_param7^2 + bias_param7^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6 + mse_param7) / 7,
    .groups = 'drop'
  )



############################################################
#### SCENARIO 5: POISSON INHOM CON xyt Z1-Z2-Z3 N = 500 ####
############################################################

inh_scen5_n500 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])}, 
                                        par = c(4.4, 2, 2, 2, -1, -1, -1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
n <- NULL
for(i in 1:length(inh_scen5_n500)) {
  n[i] <- nrow(inh_scen5_n500[[i]]$df)
}
mean(n)

true_params_scen5_n500 <- c(4.4, 2, 2, 2, -1, -1, -1)
mult_param_scen5_n500 <- seq(1, 50, by = 1)
ncube_param_scen5_n500 <- seq(1, 50, by = 1)
grid_param_scen5_n500 <- c(FALSE, TRUE)
df_param_scen5_n500 <- expand.grid(mult = mult_param_scen5_n500, ncube = ncube_param_scen5_n500, grid = grid_param_scen5_n500)

n_settings_scen5_n500 <- nrow(df_param_scen5_n500)
n_params_scen5_n500 <- length(true_params_scen5_n500)

results_list_scen5_n500 <- vector("list", length = length(inh_scen5_n500))

for (i in seq_along(inh_scen5_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen5_n500[[i]]
  
  export_list_scen5_n500 <- c("sim_process", "df_param_scen5_n500", "inh_scen5_n500", "true_params_scen5_n500", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen5_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen5_n500, function(j) {
    mult_val <- df_param_scen5_n500$mult[j]
    ncube_val <- df_param_scen5_n500$ncube[j]
    grid_val <- df_param_scen5_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + cov1 + cov2 + cov3, mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1]), est_param2 = as.numeric(est$IntCoefs[2]),
             est_param3 = as.numeric(est$IntCoefs[3]), est_param4 = as.numeric(est$IntCoefs[4]),
             est_param5 = as.numeric(est$IntCoefs[5]), est_param6 = as.numeric(est$IntCoefs[6]),
             est_param7 = as.numeric(est$IntCoefs[7])))
  })
  
  res_df <- cbind(df_param_scen5_n500, t(res_matrix))
  results_list_scen5_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen5_n500, file = "results_list_scen5_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 5 n = 500: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen5_n500), elapsed_time))
}


all_results_scen5_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen5_n500, seq_along(results_list_scen5_n500)))

summary_stats_scen5_n500 <- all_results_scen5_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    mean_param7 = mean(est_param7),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    sd_param7 = sd(est_param7),
    bias_param1 = mean_param1 - true_params_scen5_n500[1],
    bias_param2 = mean_param2 - true_params_scen5_n500[2],
    bias_param3 = mean_param3 - true_params_scen5_n500[3],
    bias_param4 = mean_param4 - true_params_scen5_n500[4],
    bias_param5 = mean_param5 - true_params_scen5_n500[5],
    bias_param6 = mean_param6 - true_params_scen5_n500[6],
    bias_param7 = mean_param7 - true_params_scen5_n500[7],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_param7 = sd_param7^2 + bias_param7^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6 + mse_param7) / 7,
    .groups = 'drop'
  )



##################################################
#### SCENARIO 6: POISSON INHOM CON m1 N = 100 ####
##################################################

set.seed(1)
g1_scen6_n100 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 2.5, nsim = 100)
set.seed(2)
g2_scen6_n100 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 3.5, nsim = 100)
set.seed(3)
g3_scen6_n100 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 4, nsim = 100)

inh_scen6_n100 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen6_n100[[i]] <- stpm(data.frame(x = c(g1_scen6_n100[[i]]$df$x, g2_scen6_n100[[i]]$df$x, g3_scen6_n100[[i]]$df$x),
                                         y = c(g1_scen6_n100[[i]]$df$y, g2_scen6_n100[[i]]$df$y, g3_scen6_n100[[i]]$df$y),
                                         t = c(g1_scen6_n100[[i]]$df$t, g2_scen6_n100[[i]]$df$t, g3_scen6_n100[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen6_n100[[i]]$df)),
                                                          rep("b", nrow(g2_scen6_n100[[i]]$df)),
                                                          rep("c", nrow(g3_scen6_n100[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen6_n100)) {
  n[i] <- nrow(inh_scen6_n100[[i]]$df)
}
mean(n)

true_params_scen6_n100 <- c(2.5, 3.5, 4)
mult_param_scen6_n100 <- seq(1, 50, by = 1)
ncube_param_scen6_n100 <- seq(1, 50, by = 1)
grid_param_scen6_n100 <- c(FALSE, TRUE)
df_param_scen6_n100 <- expand.grid(mult = mult_param_scen6_n100, ncube = ncube_param_scen6_n100, grid = grid_param_scen6_n100)

n_settings_scen6_n100 <- nrow(df_param_scen6_n100)
n_params_scen6_n100 <- length(true_params_scen6_n100)

results_list_scen6_n100 <- vector("list", length = length(inh_scen6_n100))

for (i in seq_along(inh_scen6_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen6_n100[[i]]
  
  export_list_scen6_n100 <- c("sim_process", "df_param_scen6_n100", "inh_scen6_n100", "true_params_scen6_n100", funzioni)
  clusterExport(cl, varlist = export_list_scen6_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen6_n100, function(j) {
    mult_val <- df_param_scen6_n100$mult[j]
    ncube_val <- df_param_scen6_n100$ncube[j]
    grid_val <- df_param_scen6_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[2]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[3]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen6_n100, t(res_matrix))
  results_list_scen6_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen6_n100, file = "results_list_scen6_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 6 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen6_n100), elapsed_time))
}


all_results_scen6_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen6_n100, seq_along(results_list_scen6_n100)))

summary_stats_scen6_n100 <- all_results_scen6_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    bias_param1 = mean_param1 - true_params_scen6_n100[1],
    bias_param2 = mean_param2 - true_params_scen6_n100[2],
    bias_param3 = mean_param3 - true_params_scen6_n100[3],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3) / 3,
    .groups = 'drop'
  )



##################################################
#### SCENARIO 6: POISSON INHOM CON m1 N = 250 ####
##################################################

set.seed(1)
g1_scen6_n250 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 3.5, nsim = 100)
set.seed(2)
g2_scen6_n250 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 4, nsim = 100)
set.seed(3)
g3_scen6_n250 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 5, nsim = 100)

inh_scen6_n250 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen6_n250[[i]] <- stpm(data.frame(x = c(g1_scen6_n250[[i]]$df$x, g2_scen6_n250[[i]]$df$x, g3_scen6_n250[[i]]$df$x),
                                         y = c(g1_scen6_n250[[i]]$df$y, g2_scen6_n250[[i]]$df$y, g3_scen6_n250[[i]]$df$y),
                                         t = c(g1_scen6_n250[[i]]$df$t, g2_scen6_n250[[i]]$df$t, g3_scen6_n250[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen6_n250[[i]]$df)),
                                                          rep("b", nrow(g2_scen6_n250[[i]]$df)),
                                                          rep("c", nrow(g3_scen6_n250[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen6_n250)) {
  n[i] <- nrow(inh_scen6_n250[[i]]$df)
}
mean(n)

true_params_scen6_n250 <- c(3.5, 4, 5)
mult_param_scen6_n250 <- seq(1, 50, by = 1)
ncube_param_scen6_n250 <- seq(1, 50, by = 1)
grid_param_scen6_n250 <- c(FALSE, TRUE)
df_param_scen6_n250 <- expand.grid(mult = mult_param_scen6_n250, ncube = ncube_param_scen6_n250, grid = grid_param_scen6_n250)

n_settings_scen6_n250 <- nrow(df_param_scen6_n250)
n_params_scen6_n250 <- length(true_params_scen6_n250)

results_list_scen6_n250 <- vector("list", length = length(inh_scen6_n250))

for (i in seq_along(inh_scen6_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen6_n250[[i]]
  
  export_list_scen6_n250 <- c("sim_process", "df_param_scen6_n250", "inh_scen6_n250", "true_params_scen6_n250", funzioni)
  clusterExport(cl, varlist = export_list_scen6_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen6_n250, function(j) {
    mult_val <- df_param_scen6_n250$mult[j]
    ncube_val <- df_param_scen6_n250$ncube[j]
    grid_val <- df_param_scen6_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[2]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[3]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen6_n250, t(res_matrix))
  results_list_scen6_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen6_n250, file = "results_list_scen6_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 6 n = 250: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen6_n250), elapsed_time))
}


all_results_scen6_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen6_n250, seq_along(results_list_scen6_n250)))

summary_stats_scen6_n250 <- all_results_scen6_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    bias_param1 = mean_param1 - true_params_scen6_n250[1],
    bias_param2 = mean_param2 - true_params_scen6_n250[2],
    bias_param3 = mean_param3 - true_params_scen6_n250[3],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3) / 3,
    .groups = 'drop'
  )



##################################################
#### SCENARIO 6: POISSON INHOM CON m1 N = 500 ####
##################################################

set.seed(1)
g1_scen6_n500 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 4.5, nsim = 100)
set.seed(2)
g2_scen6_n500 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 5, nsim = 100)
set.seed(3)
g3_scen6_n500 <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = 5.5, nsim = 100)

inh_scen6_n500 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen6_n500[[i]] <- stpm(data.frame(x = c(g1_scen6_n500[[i]]$df$x, g2_scen6_n500[[i]]$df$x, g3_scen6_n500[[i]]$df$x),
                                         y = c(g1_scen6_n500[[i]]$df$y, g2_scen6_n500[[i]]$df$y, g3_scen6_n500[[i]]$df$y),
                                         t = c(g1_scen6_n500[[i]]$df$t, g2_scen6_n500[[i]]$df$t, g3_scen6_n500[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen6_n500[[i]]$df)),
                                                          rep("b", nrow(g2_scen6_n500[[i]]$df)),
                                                          rep("c", nrow(g3_scen6_n500[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen6_n500)) {
  n[i] <- nrow(inh_scen6_n500[[i]]$df)
}
mean(n)

true_params_scen6_n500 <- c(4.5, 5, 5.5)
mult_param_scen6_n500 <- seq(1, 50, by = 1)
ncube_param_scen6_n500 <- seq(1, 50, by = 1)
grid_param_scen6_n500 <- c(FALSE, TRUE)
df_param_scen6_n500 <- expand.grid(mult = mult_param_scen6_n500, ncube = ncube_param_scen6_n500, grid = grid_param_scen6_n500)

n_settings_scen6_n500 <- nrow(df_param_scen6_n500)
n_params_scen6_n500 <- length(true_params_scen6_n500)

results_list_scen6_n500 <- vector("list", length = length(inh_scen6_n500))

for (i in seq_along(inh_scen6_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen6_n500[[i]]
  
  export_list_scen6_n500 <- c("sim_process", "df_param_scen6_n500", "inh_scen6_n500", "true_params_scen6_n500", funzioni)
  clusterExport(cl, varlist = export_list_scen6_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen6_n500, function(j) {
    mult_val <- df_param_scen6_n500$mult[j]
    ncube_val <- df_param_scen6_n500$ncube[j]
    grid_val <- df_param_scen6_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[2]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[3]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen6_n500, t(res_matrix))
  results_list_scen6_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen6_n500, file = "results_list_scen6_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 6 n = 500: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen6_n500), elapsed_time))
}


all_results_scen6_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen6_n500, seq_along(results_list_scen6_n500)))

summary_stats_scen6_n500 <- all_results_scen6_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    bias_param1 = mean_param1 - true_params_scen6_n500[1],
    bias_param2 = mean_param2 - true_params_scen6_n500[2],
    bias_param3 = mean_param3 - true_params_scen6_n500[3],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3) / 3,
    .groups = 'drop'
  )



########################################################
#### SCENARIO 7: POISSON INHOM CON xyt - m1 N = 100 ####
########################################################

set.seed(1)
g1_scen7_n100 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(1, 1, 1, 1), nsim = 100)
set.seed(2)
g2_scen7_n100 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(1.5, 1, 1, 1), nsim = 100)
set.seed(3)
g3_scen7_n100 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(2, 1, 1, 1), nsim = 100)

inh_scen7_n100 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen7_n100[[i]] <- stpm(data.frame(x = c(g1_scen7_n100[[i]]$df$x, g2_scen7_n100[[i]]$df$x, g3_scen7_n100[[i]]$df$x),
                                         y = c(g1_scen7_n100[[i]]$df$y, g2_scen7_n100[[i]]$df$y, g3_scen7_n100[[i]]$df$y),
                                         t = c(g1_scen7_n100[[i]]$df$t, g2_scen7_n100[[i]]$df$t, g3_scen7_n100[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen7_n100[[i]]$df)),
                                                          rep("b", nrow(g2_scen7_n100[[i]]$df)),
                                                          rep("c", nrow(g3_scen7_n100[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen7_n100)) {
  n[i] <- nrow(inh_scen7_n100[[i]]$df)
}
mean(n)

true_params_scen7_n100 <- c(1, 1.5, 2, 1, 1, 1)
mult_param_scen7_n100 <- seq(1, 50, by = 1)
ncube_param_scen7_n100 <- seq(1, 50, by = 1)
grid_param_scen7_n100 <- c(FALSE, TRUE)
df_param_scen7_n100 <- expand.grid(mult = mult_param_scen7_n100, ncube = ncube_param_scen7_n100, grid = grid_param_scen7_n100)

n_settings_scen7_n100 <- nrow(df_param_scen7_n100)
n_params_scen7_n100 <- length(true_params_scen7_n100)

results_list_scen7_n100 <- vector("list", length = length(inh_scen7_n100))

for (i in seq_along(inh_scen7_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen7_n100[[i]]
  
  export_list_scen7_n100 <- c("sim_process", "df_param_scen7_n100", "inh_scen7_n100", "true_params_scen7_n100", funzioni)
  clusterExport(cl, varlist = export_list_scen7_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen7_n100, function(j) {
    mult_val <- df_param_scen7_n100$mult[j]
    ncube_val <- df_param_scen7_n100$ncube[j]
    grid_val <- df_param_scen7_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[5]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
             est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), est_param6 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen7_n100, t(res_matrix))
  results_list_scen7_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen7_n100, file = "results_list_scen7_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 7 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen7_n100), elapsed_time))
}


all_results_scen7_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen7_n100, seq_along(results_list_scen7_n100)))

summary_stats_scen7_n100 <- all_results_scen7_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    bias_param1 = mean_param1 - true_params_scen7_n100[1],
    bias_param2 = mean_param2 - true_params_scen7_n100[2],
    bias_param3 = mean_param3 - true_params_scen7_n100[3],
    bias_param4 = mean_param4 - true_params_scen7_n100[4],
    bias_param5 = mean_param5 - true_params_scen7_n100[5],
    bias_param6 = mean_param6 - true_params_scen7_n100[6],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6) / 6,
    .groups = 'drop'
  )



########################################################
#### SCENARIO 7: POISSON INHOM CON xyt - m1 N = 250 ####
########################################################

set.seed(1)
g1_scen7_n250 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(2, 1, 1, 1), nsim = 100)
set.seed(2)
g2_scen7_n250 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(2.5, 1, 1, 1), nsim = 100)
set.seed(3)
g3_scen7_n250 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(3, 1, 1, 1), nsim = 100)

inh_scen7_n250 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen7_n250[[i]] <- stpm(data.frame(x = c(g1_scen7_n250[[i]]$df$x, g2_scen7_n250[[i]]$df$x, g3_scen7_n250[[i]]$df$x),
                                         y = c(g1_scen7_n250[[i]]$df$y, g2_scen7_n250[[i]]$df$y, g3_scen7_n250[[i]]$df$y),
                                         t = c(g1_scen7_n250[[i]]$df$t, g2_scen7_n250[[i]]$df$t, g3_scen7_n250[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen7_n250[[i]]$df)),
                                                          rep("b", nrow(g2_scen7_n250[[i]]$df)),
                                                          rep("c", nrow(g3_scen7_n250[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen7_n250)) {
  n[i] <- nrow(inh_scen7_n250[[i]]$df)
}
mean(n)

true_params_scen7_n250 <- c(2, 2.5, 3, 1, 1, 1)
mult_param_scen7_n250 <- seq(1, 50, by = 1)
ncube_param_scen7_n250 <- seq(1, 50, by = 1)
grid_param_scen7_n250 <- c(FALSE, TRUE)
df_param_scen7_n250 <- expand.grid(mult = mult_param_scen7_n250, ncube = ncube_param_scen7_n250, grid = grid_param_scen7_n250)

n_settings_scen7_n250 <- nrow(df_param_scen7_n250)
n_params_scen7_n250 <- length(true_params_scen7_n250)

results_list_scen7_n250 <- vector("list", length = length(inh_scen7_n250))

for (i in seq_along(inh_scen7_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen7_n250[[i]]
  
  export_list_scen7_n250 <- c("sim_process", "df_param_scen7_n250", "inh_scen7_n250", "true_params_scen7_n250", funzioni)
  clusterExport(cl, varlist = export_list_scen7_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen7_n250, function(j) {
    mult_val <- df_param_scen7_n250$mult[j]
    ncube_val <- df_param_scen7_n250$ncube[j]
    grid_val <- df_param_scen7_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[5]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
             est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), est_param6 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen7_n250, t(res_matrix))
  results_list_scen7_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen7_n250, file = "results_list_scen7_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 7 n = 250: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen7_n250), elapsed_time))
}


all_results_scen7_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen7_n250, seq_along(results_list_scen7_n250)))

summary_stats_scen7_n250 <- all_results_scen7_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    bias_param1 = mean_param1 - true_params_scen7_n250[1],
    bias_param2 = mean_param2 - true_params_scen7_n250[2],
    bias_param3 = mean_param3 - true_params_scen7_n250[3],
    bias_param4 = mean_param4 - true_params_scen7_n250[4],
    bias_param5 = mean_param5 - true_params_scen7_n250[5],
    bias_param6 = mean_param6 - true_params_scen7_n250[6],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6) / 6,
    .groups = 'drop'
  )



########################################################
#### SCENARIO 7: POISSON INHOM CON xyt - m1 N = 500 ####
########################################################

set.seed(1)
g1_scen7_n500 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(2.5, 1, 1, 1), nsim = 100)
set.seed(2)
g2_scen7_n500 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(3, 1, 1, 1), nsim = 100)
set.seed(3)
g3_scen7_n500 <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(4, 1, 1, 1), nsim = 100)

inh_scen7_n500 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen7_n500[[i]] <- stpm(data.frame(x = c(g1_scen7_n500[[i]]$df$x, g2_scen7_n500[[i]]$df$x, g3_scen7_n500[[i]]$df$x),
                                         y = c(g1_scen7_n500[[i]]$df$y, g2_scen7_n500[[i]]$df$y, g3_scen7_n500[[i]]$df$y),
                                         t = c(g1_scen7_n500[[i]]$df$t, g2_scen7_n500[[i]]$df$t, g3_scen7_n500[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen7_n500[[i]]$df)),
                                                          rep("b", nrow(g2_scen7_n500[[i]]$df)),
                                                          rep("c", nrow(g3_scen7_n500[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen7_n500)) {
  n[i] <- nrow(inh_scen7_n500[[i]]$df)
}
mean(n)

true_params_scen7_n500 <- c(2.5, 3, 4, 1, 1, 1)
mult_param_scen7_n500 <- seq(1, 50, by = 1)
ncube_param_scen7_n500 <- seq(1, 50, by = 1)
grid_param_scen7_n500 <- c(FALSE, TRUE)
df_param_scen7_n500 <- expand.grid(mult = mult_param_scen7_n500, ncube = ncube_param_scen7_n500, grid = grid_param_scen7_n500)

n_settings_scen7_n500 <- nrow(df_param_scen7_n500)
n_params_scen7_n500 <- length(true_params_scen7_n500)

results_list_scen7_n500 <- vector("list", length = length(inh_scen7_n500))

for (i in seq_along(inh_scen7_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen7_n500[[i]]
  
  export_list_scen7_n500 <- c("sim_process", "df_param_scen7_n500", "inh_scen7_n500", "true_params_scen7_n500", funzioni)
  clusterExport(cl, varlist = export_list_scen7_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen7_n500, function(j) {
    mult_val <- df_param_scen7_n500$mult[j]
    ncube_val <- df_param_scen7_n500$ncube[j]
    grid_val <- df_param_scen7_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[5]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
             est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), est_param6 = as.numeric(est$IntCoefs[4])))
  })
  
  res_df <- cbind(df_param_scen7_n500, t(res_matrix))
  results_list_scen7_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen7_n500, file = "results_list_scen7_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 7 n = 500: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen7_n500), elapsed_time))
}


all_results_scen7_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen7_n500, seq_along(results_list_scen7_n500)))

summary_stats_scen7_n500 <- all_results_scen7_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    bias_param1 = mean_param1 - true_params_scen7_n500[1],
    bias_param2 = mean_param2 - true_params_scen7_n500[2],
    bias_param3 = mean_param3 - true_params_scen7_n500[3],
    bias_param4 = mean_param4 - true_params_scen7_n500[4],
    bias_param5 = mean_param5 - true_params_scen7_n500[5],
    bias_param6 = mean_param6 - true_params_scen7_n500[6],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6) / 6,
    .groups = 'drop'
  )



#############################################################
#### SCENARIO 8: POISSON INHOM CON xyt - m1 - Z1 N = 100 ####
#############################################################

g1_scen8_n100 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(1, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 1)
g2_scen8_n100 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(1.3, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
g3_scen8_n100 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(1.5, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 3)

inh_scen8_n100 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen8_n100[[i]] <- stpm(data.frame(x = c(g1_scen8_n100[[i]]$df$x, g2_scen8_n100[[i]]$df$x, g3_scen8_n100[[i]]$df$x),
                                         y = c(g1_scen8_n100[[i]]$df$y, g2_scen8_n100[[i]]$df$y, g3_scen8_n100[[i]]$df$y),
                                         t = c(g1_scen8_n100[[i]]$df$t, g2_scen8_n100[[i]]$df$t, g3_scen8_n100[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen8_n100[[i]]$df)),
                                                          rep("b", nrow(g2_scen8_n100[[i]]$df)),
                                                          rep("c", nrow(g3_scen8_n100[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen8_n100)) {
  n[i] <- nrow(inh_scen8_n100[[i]]$df)
}
mean(n)

true_params_scen8_n100 <- c(1, 1.3, 1.5, 1, 1, 1, 1)
mult_param_scen8_n100 <- seq(1, 50, by = 1)
ncube_param_scen8_n100 <- seq(1, 50, by = 1)
grid_param_scen8_n100 <- c(FALSE, TRUE)
df_param_scen8_n100 <- expand.grid(mult = mult_param_scen8_n100, ncube = ncube_param_scen8_n100, grid = grid_param_scen8_n100)

n_settings_scen8_n100 <- nrow(df_param_scen8_n100)
n_params_scen8_n100 <- length(true_params_scen8_n100)

results_list_scen8_n100 <- vector("list", length = length(inh_scen8_n100))

for (i in seq_along(inh_scen8_n100)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen8_n100[[i]]
  
  export_list_scen8_n100 <- c("sim_process", "df_param_scen8_n100", "inh_scen8_n100", "true_params_scen8_n100", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen8_n100)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen8_n100, function(j) {
    mult_val <- df_param_scen8_n100$mult[j]
    ncube_val <- df_param_scen8_n100$ncube[j]
    grid_val <- df_param_scen8_n100$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[8]),
             est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), 
             est_param6 = as.numeric(est$IntCoefs[4]), est_param7 = as.numeric(est$IntCoefs[5])))
  })
  
  res_df <- cbind(df_param_scen8_n100, t(res_matrix))
  results_list_scen8_n100[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen8_n100, file = "results_list_scen8_n100.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 8 n = 100: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen8_n100), elapsed_time))
}


all_results_scen8_n100 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen8_n100, seq_along(results_list_scen8_n100)))

summary_stats_scen8_n100 <- all_results_scen8_n100 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    mean_param7 = mean(est_param7),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    sd_param7 = sd(est_param7),
    bias_param1 = mean_param1 - true_params_scen8_n100[1],
    bias_param2 = mean_param2 - true_params_scen8_n100[2],
    bias_param3 = mean_param3 - true_params_scen8_n100[3],
    bias_param4 = mean_param4 - true_params_scen8_n100[4],
    bias_param5 = mean_param5 - true_params_scen8_n100[5],
    bias_param6 = mean_param6 - true_params_scen8_n100[6],
    bias_param7 = mean_param7 - true_params_scen8_n100[7],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_param7 = sd_param7^2 + bias_param7^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6 + mse_param7) / 7,
    .groups = 'drop'
  )



#############################################################
#### SCENARIO 8: POISSON INHOM CON xyt - m1 - Z1 N = 250 ####
#############################################################

g1_scen8_n250 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(1.8, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 1)
g2_scen8_n250 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(2.2, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
g3_scen8_n250 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(2.5, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 3)

inh_scen8_n250 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen8_n250[[i]] <- stpm(data.frame(x = c(g1_scen8_n250[[i]]$df$x, g2_scen8_n250[[i]]$df$x, g3_scen8_n250[[i]]$df$x),
                                         y = c(g1_scen8_n250[[i]]$df$y, g2_scen8_n250[[i]]$df$y, g3_scen8_n250[[i]]$df$y),
                                         t = c(g1_scen8_n250[[i]]$df$t, g2_scen8_n250[[i]]$df$t, g3_scen8_n250[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen8_n250[[i]]$df)),
                                                          rep("b", nrow(g2_scen8_n250[[i]]$df)),
                                                          rep("c", nrow(g3_scen8_n250[[i]]$df))))))
}

n <- NULL
for(i in 1:length(inh_scen8_n250)) {
  n[i] <- nrow(inh_scen8_n250[[i]]$df)
}
mean(n)

true_params_scen8_n250 <- c(1.8, 2.2, 2.5, 1, 1, 1, 1)
mult_param_scen8_n250 <- seq(1, 50, by = 1)
ncube_param_scen8_n250 <- seq(1, 50, by = 1)
grid_param_scen8_n250 <- c(FALSE, TRUE)
df_param_scen8_n250 <- expand.grid(mult = mult_param_scen8_n250, ncube = ncube_param_scen8_n250, grid = grid_param_scen8_n250)

n_settings_scen8_n250 <- nrow(df_param_scen8_n250)
n_params_scen8_n250 <- length(true_params_scen8_n250)

results_list_scen8_n250 <- vector("list", length = length(inh_scen8_n250))

for (i in seq_along(inh_scen8_n250)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen8_n250[[i]]
  
  export_list_scen8_n250 <- c("sim_process", "df_param_scen8_n250", "inh_scen8_n250", "true_params_scen8_n250", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen8_n250)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen8_n250, function(j) {
    mult_val <- df_param_scen8_n250$mult[j]
    ncube_val <- df_param_scen8_n250$ncube[j]
    grid_val <- df_param_scen8_n250$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[8]),
             est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), 
             est_param6 = as.numeric(est$IntCoefs[4]), est_param7 = as.numeric(est$IntCoefs[5])))
  })
  
  res_df <- cbind(df_param_scen8_n250, t(res_matrix))
  results_list_scen8_n250[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen8_n250, file = "results_list_scen8_n250.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 8 n = 250: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen8_n250), elapsed_time))
}


all_results_scen8_n250 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen8_n250, seq_along(results_list_scen8_n250)))

summary_stats_scen8_n250 <- all_results_scen8_n250 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    mean_param7 = mean(est_param7),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    sd_param7 = sd(est_param7),
    bias_param1 = mean_param1 - true_params_scen8_n250[1],
    bias_param2 = mean_param2 - true_params_scen8_n250[2],
    bias_param3 = mean_param3 - true_params_scen8_n250[3],
    bias_param4 = mean_param4 - true_params_scen8_n250[4],
    bias_param5 = mean_param5 - true_params_scen8_n250[5],
    bias_param6 = mean_param6 - true_params_scen8_n250[6],
    bias_param7 = mean_param7 - true_params_scen8_n250[7],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_param7 = sd_param7^2 + bias_param7^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6 + mse_param7) / 7,
    .groups = 'drop'
  )



#############################################################
#### SCENARIO 8: POISSON INHOM CON xyt - m1 - Z1 N = 500 ####
#############################################################

g1_scen8_n500 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(2.5, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 1)
g2_scen8_n500 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(3, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 2)
g3_scen8_n500 <- rstpp_cov_generalised(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(3.2, 1, 1, 1, 1), covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3), nsim = 100, seed = 3)

inh_scen8_n500 <- setNames(replicate(100, data.frame()), paste0("X", 1:100))

for(i in 1:100) {
  inh_scen8_n500[[i]] <- stpm(data.frame(x = c(g1_scen8_n500[[i]]$df$x, g2_scen8_n500[[i]]$df$x, g3_scen8_n500[[i]]$df$x),
                                         y = c(g1_scen8_n500[[i]]$df$y, g2_scen8_n500[[i]]$df$y, g3_scen8_n500[[i]]$df$y),
                                         t = c(g1_scen8_n500[[i]]$df$t, g2_scen8_n500[[i]]$df$t, g3_scen8_n500[[i]]$df$t),
                                         m1 = as.factor(c(rep("a", nrow(g1_scen8_n500[[i]]$df)),
                                                          rep("b", nrow(g2_scen8_n500[[i]]$df)),
                                                          rep("c", nrow(g3_scen8_n500[[i]]$df))))))
}


inh_scen8_n500[[1]]
est <- stppm_prova_local(X = inh_scen8_n500[[1]], formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                         marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)


n <- NULL
for(i in 1:length(inh_scen8_n500)) {
  n[i] <- nrow(inh_scen8_n500[[i]]$df)
}
mean(n)

true_params_scen8_n500 <- c(2.5, 3, 3.2, 1, 1, 1, 1)
mult_param_scen8_n500 <- seq(1, 50, by = 1)
ncube_param_scen8_n500 <- seq(1, 50, by = 1)
grid_param_scen8_n500 <- c(FALSE, TRUE)
df_param_scen8_n500 <- expand.grid(mult = mult_param_scen8_n500, ncube = ncube_param_scen8_n500, grid = grid_param_scen8_n500)

n_settings_scen8_n500 <- nrow(df_param_scen8_n500)
n_params_scen8_n500 <- length(true_params_scen8_n500)

results_list_scen8_n500 <- vector("list", length = length(inh_scen8_n500))

for (i in seq_along(inh_scen8_n500)) {
  
  start_time <- Sys.time()
  sim_process <- inh_scen8_n500[[i]]
  
  export_list_scen8_n500 <- c("sim_process", "df_param_scen8_n500", "inh_scen8_n500", "true_params_scen8_n500", "df_covs", funzioni)
  clusterExport(cl, varlist = export_list_scen8_n500)
  
  res_matrix <- parSapply(cl, 1:n_settings_scen8_n500, function(j) {
    mult_val <- df_param_scen8_n500$mult[j]
    ncube_val <- df_param_scen8_n500$ncube[j]
    grid_val <- df_param_scen8_n500$grid[j]
    
    est <- stppm_prova_local(X = sim_process, formula = ~ x + y + t + cov1 + s(m1, bs = "re"), mult = mult_val, ncube = ncube_val, grid = grid_val, 
                             marked = TRUE, spatial.cov = TRUE, covs = df_covs, seed = 2)
    
    return(c(est_param1 = as.numeric(est$IntCoefs[1] + est$IntCoefs[6]), 
             est_param2 = as.numeric(est$IntCoefs[1] + est$IntCoefs[7]),
             est_param3 = as.numeric(est$IntCoefs[1] + est$IntCoefs[8]),
             est_param4 = as.numeric(est$IntCoefs[2]), est_param5 = as.numeric(est$IntCoefs[3]), 
             est_param6 = as.numeric(est$IntCoefs[4]), est_param7 = as.numeric(est$IntCoefs[5])))
  })
  
  res_df <- cbind(df_param_scen8_n500, t(res_matrix))
  results_list_scen8_n500[[i]] <- as.data.frame(res_df)
  
  save(results_list_scen8_n500, file = "results_list_scen8_n500.RData")
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf(">> Scen 8 n = 500: Completata simulazione %d su %d (%.2f secondi)\n", i, length(inh_scen8_n500), elapsed_time))
}


all_results_scen8_n500 <- do.call(rbind, Map(function(df, i) {
  df$sim <- i
  return(df)
}, results_list_scen8_n500, seq_along(results_list_scen8_n500)))

summary_stats_scen8_n500 <- all_results_scen8_n500 %>%
  group_by(mult, ncube, grid) %>%
  summarise(
    mean_param1 = mean(est_param1),
    mean_param2 = mean(est_param2),
    mean_param3 = mean(est_param3),
    mean_param4 = mean(est_param4),
    mean_param5 = mean(est_param5),
    mean_param6 = mean(est_param6),
    mean_param7 = mean(est_param7),
    sd_param1 = sd(est_param1),
    sd_param2 = sd(est_param2),
    sd_param3 = sd(est_param3),
    sd_param4 = sd(est_param4),
    sd_param5 = sd(est_param5),
    sd_param6 = sd(est_param6),
    sd_param7 = sd(est_param7),
    bias_param1 = mean_param1 - true_params_scen8_n500[1],
    bias_param2 = mean_param2 - true_params_scen8_n500[2],
    bias_param3 = mean_param3 - true_params_scen8_n500[3],
    bias_param4 = mean_param4 - true_params_scen8_n500[4],
    bias_param5 = mean_param5 - true_params_scen8_n500[5],
    bias_param6 = mean_param6 - true_params_scen8_n500[6],
    bias_param7 = mean_param7 - true_params_scen8_n500[7],
    mse_param1 = sd_param1^2 + bias_param1^2,
    mse_param2 = sd_param2^2 + bias_param2^2,
    mse_param3 = sd_param3^2 + bias_param3^2,
    mse_param4 = sd_param4^2 + bias_param4^2,
    mse_param5 = sd_param5^2 + bias_param5^2,
    mse_param6 = sd_param6^2 + bias_param6^2,
    mse_param7 = sd_param7^2 + bias_param7^2,
    mse_combined = (mse_param1 + mse_param2 + mse_param3 + mse_param4 + mse_param5 + mse_param6 + mse_param7) / 7,
    .groups = 'drop'
  )
