

######################################
#### CODES SIMULATION STUDY - GET ####
######################################

set_contour_class_na <- function(df, combinations_to_na) {
  # cicla sulle combinazioni da sostituire
  for (i in seq_len(nrow(combinations_to_na))) {
    m <- combinations_to_na$mult[i]
    n <- combinations_to_na$ncube[i]
    
    df$contour_class[df$mult == m & df$ncube == n] <- NA
  }
  return(df)
}


#### DO NOT RUN: COMPUTATIONALLY VERY EXPENSIVE ####



########################
### SCEN 1 - N = 100 ###
########################

true_params_scen1_n100
df_summary_stats_gridF_scen1_n100 <- df_summary_stats_scen1_n100[df_summary_stats_scen1_n100$grid == FALSE,]
df_summary_stats_gridT_scen1_n100 <- df_summary_stats_scen1_n100[df_summary_stats_scen1_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen1_n100 <- df_summary_stats_gridF_scen1_n100
df_summary_stats_pvalue_gridF_scen1_n100$pvalue <- NA

min_mse_gridF_scen1_n100 <- min(df_summary_stats_gridF_scen1_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen1_n100)
count_significant <- 0
prop_explored <- 0

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen1_n100[[id]]

ncores <- 80
cl <- makeCluster(ncores)
registerDoParallel(cl)

funzioni_da_esportare <- c("interp3D_lapply_mod", "stp", "Khat_spatial.3D", 
                           "interp3D_point_mod", "stikfunction.spatial.3D", 
                           "edge.Trans.4D", "sbox.4D", "rstpp_cov_generalised_v2", "interp3D_prova3")

clusterExport(cl, varlist = funzioni_da_esportare)

clusterEvalQ(cl, {
  library(spatstat)
  library(stpp) 
  library(stopp)
  library(mgcv)
  library(GET)
})

while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen1_n100 <- min_mse_gridF_scen1_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen1_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen1_n100 & is.na(pvalue))
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x)},
                         par = c(comb$mean_param1, comb$mean_param2), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen1_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen1_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen1_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen1_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen1_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen1_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen1_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen1_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}

stopCluster(cl)

df_summary_stats_pvalue_gridF_scen1_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen1_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

library(ggplot2)
library(grid)

library(ggplot2)
library(grid)

p_gridF_scen1_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen1_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen1_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen1_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen1_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen1_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen1_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen1_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 14 | nc = 4"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen1_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen1_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen1_n100)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))
comb_subset[which.max(comb_subset$mse_combined),]

save(df_summary_stats_pvalue_gridF_scen1_n100, file = "pvalue_plot_gridF_scen1_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen1_n100 <- df_summary_stats_gridT_scen1_n100
df_summary_stats_pvalue_gridT_scen1_n100$pvalue <- NA

min_mse_gridT_scen1_n100 <- min(df_summary_stats_gridT_scen1_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen1_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen1_n100[[id]]

ncores <- 64
cl <- makeCluster(ncores)
registerDoParallel(cl)

funzioni_da_esportare <- c("interp3D_lapply_mod", "stp", "Khat_spatial.3D", 
                           "interp3D_point_mod", "stikfunction.spatial.3D", 
                           "edge.Trans.4D", "sbox.4D")

clusterExport(cl, varlist = funzioni_da_esportare)

clusterEvalQ(cl, {
  library(spatstat)
  library(stpp) 
  library(stopp)
  library(mgcv)
  library(GET)
})

while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen1_n100 <- min_mse_gridT_scen1_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen1_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen1_n100 & is.na(pvalue))
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x)},
                         par = c(comb$mean_param1, comb$mean_param2), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen1_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen1_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen1_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen1_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen1_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen1_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen1_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen1_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}

stopCluster(cl)

df_summary_stats_pvalue_gridT_scen1_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen1_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen1_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen1_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen1_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen1_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen1_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen1_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen1_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen1_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 9 | nc = 7"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen1_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen1_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen1_n100)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen1_n100, file = "pvalue_plot_gridT_scen1_n100.RData")





########################
### SCEN 1 - N = 250 ###
########################

true_params_scen1_n250
df_summary_stats_gridF_scen1_n250 <- df_summary_stats_scen1_n250[df_summary_stats_scen1_n250$grid == FALSE,]
df_summary_stats_gridT_scen1_n250 <- df_summary_stats_scen1_n250[df_summary_stats_scen1_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen1_n250 <- df_summary_stats_gridF_scen1_n250
df_summary_stats_pvalue_gridF_scen1_n250$pvalue <- NA

min_mse_gridF_scen1_n250 <- min(df_summary_stats_gridF_scen1_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen1_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen1_n250[[id]]

ncores <- 64
cl <- makeCluster(ncores)
registerDoParallel(cl)

funzioni_da_esportare <- c("interp3D_lapply_mod", "stp", "Khat_spatial.3D", 
                           "interp3D_point_mod", "stikfunction.spatial.3D", 
                           "edge.Trans.4D", "sbox.4D")

clusterExport(cl, varlist = funzioni_da_esportare)

clusterEvalQ(cl, {
  library(spatstat)
  library(stpp) 
  library(stopp)
  library(mgcv)
  library(GET)
})

while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen1_n250 <- min_mse_gridF_scen1_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen1_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen1_n250 & is.na(pvalue))
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x)},
                         par = c(comb$mean_param1, comb$mean_param2), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen1_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen1_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen1_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen1_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen1_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen1_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen1_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen1_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}

stopCluster(cl)

df_summary_stats_pvalue_gridF_scen1_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen1_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen1_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen1_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen1_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen1_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen1_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen1_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen1_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen1_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 10 | nc = 6"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen1_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen1_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen1_n250)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen1_n250, file = "pvalue_plot_gridF_scen1_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen1_n250 <- df_summary_stats_gridT_scen1_n250
df_summary_stats_pvalue_gridT_scen1_n250$pvalue <- NA

min_mse_gridT_scen1_n250 <- min(df_summary_stats_gridT_scen1_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen1_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen1_n250[[id]]

ncores <- 64
cl <- makeCluster(ncores)
registerDoParallel(cl)

funzioni_da_esportare <- c("interp3D_lapply_mod", "stp", "Khat_spatial.3D", 
                           "interp3D_point_mod", "stikfunction.spatial.3D", 
                           "edge.Trans.4D", "sbox.4D")

clusterExport(cl, varlist = funzioni_da_esportare)

clusterEvalQ(cl, {
  library(spatstat)
  library(stpp) 
  library(stopp)
  library(mgcv)
  library(GET)
})

while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen1_n250 <- min_mse_gridT_scen1_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen1_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen1_n250 & is.na(pvalue))
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x)},
                         par = c(comb$mean_param1, comb$mean_param2), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen1_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen1_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen1_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen1_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen1_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen1_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen1_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen1_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}

stopCluster(cl)

df_summary_stats_pvalue_gridT_scen1_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen1_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen1_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen1_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen1_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen1_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen1_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen1_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen1_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen1_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 38 | nc = 17"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen1_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen1_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen1_n250)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen1_n250, file = "pvalue_plot_gridT_scen1_n250.RData")





########################
### SCEN 1 - N = 500 ###
########################

true_params_scen1_n500
df_summary_stats_gridF_scen1_n500 <- df_summary_stats_scen1_n500[df_summary_stats_scen1_n500$grid == FALSE,]
df_summary_stats_gridT_scen1_n500 <- df_summary_stats_scen1_n500[df_summary_stats_scen1_n500$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen1_n500 <- df_summary_stats_gridF_scen1_n500
df_summary_stats_pvalue_gridF_scen1_n500$pvalue <- NA

min_mse_gridF_scen1_n500 <- min(df_summary_stats_gridF_scen1_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen1_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen1_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen1_n500 <- min_mse_gridF_scen1_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen1_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen1_n500 & is.na(pvalue))
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x)},
                         par = c(comb$mean_param1, comb$mean_param2), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen1_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen1_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen1_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen1_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen1_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen1_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen1_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen1_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen1_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen1_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen1_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen1_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen1_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen1_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen1_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen1_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen1_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen1_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 5 | nc = 6"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen1_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen1_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen1_n500)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen1_n500, file = "pvalue_plot_gridF_scen1_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen1_n500 <- df_summary_stats_gridT_scen1_n500
df_summary_stats_pvalue_gridT_scen1_n500$pvalue <- NA

min_mse_gridT_scen1_n500 <- min(df_summary_stats_gridT_scen1_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen1_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen1_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen1_n500 <- min_mse_gridT_scen1_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen1_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen1_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen1_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen1_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    comb <- comb_subset[i, ]
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x)},
                         par = c(comb$mean_param1, comb$mean_param2), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen1_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen1_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen1_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen1_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen1_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen1_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  
  # Incremento soglia per il prossimo giro
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen1_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen1_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen1_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen1_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen1_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen1_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen1_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen1_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen1_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen1_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen1_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 40 | nc = 23"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen1_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen1_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen1_n500)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen1_n500, file = "pvalue_plot_gridT_scen1_n500.RData")





########################
### SCEN 2 - N = 100 ###
########################

true_params_scen2_n100
df_summary_stats_gridF_scen2_n100 <- df_summary_stats_scen2_n100[df_summary_stats_scen2_n100$grid == FALSE,]
df_summary_stats_gridT_scen2_n100 <- df_summary_stats_scen2_n100[df_summary_stats_scen2_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen2_n100 <- df_summary_stats_gridF_scen2_n100
df_summary_stats_pvalue_gridF_scen2_n100$pvalue <- NA

min_mse_gridF_scen2_n100 <- min(df_summary_stats_gridF_scen2_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen2_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen2_n100[[id+2]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen2_n100 <- min_mse_gridF_scen2_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen2_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen2_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen2_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen2_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x +
                            comb$mean_param3 * false_ppp$df$y + comb$mean_param4 * false_ppp$df$t)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x  + a[3] * y  + a[4] * t)},
                         par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x +
                        comb$mean_param3 * bestmod_ppp[[j]]$df$y + comb$mean_param4 * bestmod_ppp[[j]]$df$t)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen2_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen2_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen2_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen2_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen2_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen2_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen2_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen2_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen2_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen2_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen2_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen2_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen2_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen2_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen2_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen2_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen2_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen2_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 48 | nc = 8"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen2_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen2_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen2_n100)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen2_n100, file = "pvalue_plot_gridF_scen2_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen2_n100 <- df_summary_stats_gridT_scen2_n100
df_summary_stats_pvalue_gridT_scen2_n100$pvalue <- NA

min_mse_gridT_scen2_n100 <- min(df_summary_stats_gridT_scen2_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen2_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen2_n100[[id+2]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen2_n100 <- min_mse_gridT_scen2_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen2_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen2_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen2_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen2_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x + 
                           comb$mean_param3 * true_ppp$df$y + comb$mean_param4 * true_ppp$df$t)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x + a[3] * y + a[4] * t)},
                         par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x +
                        comb$mean_param3 * bestmod_ppp[[j]]$df$y + comb$mean_param4 * bestmod_ppp[[j]]$df$t)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen2_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen2_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen2_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen2_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen2_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen2_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen2_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen2_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen2_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen2_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen2_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen2_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen2_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen2_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen2_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen2_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen2_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen2_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 50 | nc = 15"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen2_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen2_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen2_n100)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen2_n100, file = "pvalue_plot_gridT_scen2_n100.RData")





########################
### SCEN 2 - N = 250 ###
########################

true_params_scen2_n250
df_summary_stats_gridF_scen2_n250 <- df_summary_stats_scen2_n250[df_summary_stats_scen2_n250$grid == FALSE,]
df_summary_stats_gridT_scen2_n250 <- df_summary_stats_scen2_n250[df_summary_stats_scen2_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen2_n250 <- df_summary_stats_gridF_scen2_n250
df_summary_stats_pvalue_gridF_scen2_n250$pvalue <- NA

min_mse_gridF_scen2_n250 <- min(df_summary_stats_gridF_scen2_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen2_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen2_n250[[id+1]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen2_n250 <- min_mse_gridF_scen2_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen2_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen2_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen2_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen2_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x +
                            comb$mean_param3 * false_ppp$df$y + comb$mean_param4 * false_ppp$df$t)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x  + a[3] * y  + a[4] * t)},
                         par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x +
                        comb$mean_param3 * bestmod_ppp[[j]]$df$y + comb$mean_param4 * bestmod_ppp[[j]]$df$t)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen2_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen2_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen2_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen2_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen2_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen2_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen2_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen2_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen2_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen2_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen2_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen2_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen2_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen2_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen2_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen2_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen2_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen2_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 42 | nc = 11"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen2_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen2_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen2_n250)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen2_n250, file = "pvalue_plot_gridF_scen2_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen2_n250 <- df_summary_stats_gridT_scen2_n250
df_summary_stats_pvalue_gridT_scen2_n250$pvalue <- NA

min_mse_gridT_scen2_n250 <- min(df_summary_stats_gridT_scen2_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen2_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen2_n250[[id+1]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen2_n250 <- min_mse_gridT_scen2_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen2_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen2_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen2_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen2_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x + 
                           comb$mean_param3 * true_ppp$df$y + comb$mean_param4 * true_ppp$df$t)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x + a[3] * y + a[4] * t)},
                         par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x +
                        comb$mean_param3 * bestmod_ppp[[j]]$df$y + comb$mean_param4 * bestmod_ppp[[j]]$df$t)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen2_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen2_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen2_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen2_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen2_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen2_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen2_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen2_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen2_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen2_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen2_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen2_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen2_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen2_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen2_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen2_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen2_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen2_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 23 | nc = 15"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen2_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen2_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen2_n250)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen2_n250, file = "pvalue_plot_gridT_scen2_n250.RData")





########################
### SCEN 2 - N = 500 ###
########################

true_params_scen2_n500
df_summary_stats_gridF_scen2_n500 <- df_summary_stats_scen2_n500[df_summary_stats_scen2_n500$grid == FALSE,]
df_summary_stats_gridT_scen2_n500 <- df_summary_stats_scen2_n500[df_summary_stats_scen2_n500$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen2_n500 <- df_summary_stats_gridF_scen2_n500
df_summary_stats_pvalue_gridF_scen2_n500$pvalue <- NA

min_mse_gridF_scen2_n500 <- min(df_summary_stats_gridF_scen2_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen2_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen2_n500[[id+2]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen2_n500 <- min_mse_gridF_scen2_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen2_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen2_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen2_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen2_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x +
                            comb$mean_param3 * false_ppp$df$y + comb$mean_param4 * false_ppp$df$t)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x  + a[3] * y  + a[4] * t)},
                         par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x +
                        comb$mean_param3 * bestmod_ppp[[j]]$df$y + comb$mean_param4 * bestmod_ppp[[j]]$df$t)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen2_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen2_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen2_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen2_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen2_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen2_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen2_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen2_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen2_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen2_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen2_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen2_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen2_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen2_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen2_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen2_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen2_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen2_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 41 | nc = 13"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen2_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen2_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen2_n500)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen2_n500, file = "pvalue_plot_gridF_scen2_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen2_n500 <- df_summary_stats_gridT_scen2_n500
df_summary_stats_pvalue_gridT_scen2_n500$pvalue <- NA

min_mse_gridT_scen2_n500 <- min(df_summary_stats_gridT_scen2_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen2_n500)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen2_n500[[id+2]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen2_n500 <- min_mse_gridT_scen2_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen2_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen2_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen2_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen2_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x + 
                           comb$mean_param3 * true_ppp$df$y + comb$mean_param4 * true_ppp$df$t)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    
    set.seed(2)
    bestmod_ppp <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2] * x + a[3] * y + a[4] * t)},
                         par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100)
    
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = 100))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$df$x +
                        comb$mean_param3 * bestmod_ppp[[j]]$df$y + comb$mean_param4 * bestmod_ppp[[j]]$df$t)
      Khat_sim_j <- Khat_spatial.3D(bestmod_ppp[[j]], lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen2_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen2_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen2_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen2_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen2_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen2_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen2_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen2_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen2_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen2_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen2_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen2_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen2_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen2_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen2_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen2_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen2_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen2_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 39 | nc = 23"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen2_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen2_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen2_n500)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen2_n500, file = "pvalue_plot_gridT_scen2_n500.RData")





########################
### SCEN 3 - N = 100 ###
########################

true_params_scen3_n100
df_summary_stats_gridF_scen3_n100 <- df_summary_stats_scen3_n100[df_summary_stats_scen3_n100$grid == FALSE,]
df_summary_stats_gridT_scen3_n100 <- df_summary_stats_scen3_n100[df_summary_stats_scen3_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen3_n100 <- df_summary_stats_gridF_scen3_n100
df_summary_stats_pvalue_gridF_scen3_n100$pvalue <- NA

min_mse_gridF_scen3_n100 <- min(df_summary_stats_gridF_scen3_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen3_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen3_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen3_n100 <- min_mse_gridF_scen3_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen3_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen3_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen3_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen3_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * covariate_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                            par = c(comb$mean_param1, comb$mean_param2), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen3_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen3_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen3_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen3_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen3_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen3_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen3_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen3_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen3_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen3_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen3_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen3_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen3_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen3_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen3_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen3_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen3_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen3_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 6 | nc = 4"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen3_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen3_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen3_n100)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen3_n100, file = "pvalue_plot_gridF_scen3_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen3_n100 <- df_summary_stats_gridT_scen3_n100
df_summary_stats_pvalue_gridT_scen3_n100$pvalue <- NA

min_mse_gridT_scen3_n100 <- min(df_summary_stats_gridT_scen3_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen3_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen3_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen3_n100 <- min_mse_gridT_scen3_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen3_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen3_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen3_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen3_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * covariate_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                            par = c(comb$mean_param1, comb$mean_param2), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen3_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen3_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen3_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen3_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen3_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen3_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen3_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen3_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen3_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen3_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen3_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen3_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen3_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen3_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen3_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen3_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen3_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen3_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 17 | nc = 7"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen3_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen3_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen3_n100)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridT_scen3_n100, file = "pvalue_plot_gridT_scen3_n100.RData")





########################
### SCEN 3 - N = 250 ###
########################

true_params_scen3_n250
df_summary_stats_gridF_scen3_n250 <- df_summary_stats_scen3_n250[df_summary_stats_scen3_n250$grid == FALSE,]
df_summary_stats_gridT_scen3_n250 <- df_summary_stats_scen3_n250[df_summary_stats_scen3_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen3_n250 <- df_summary_stats_gridF_scen3_n250
df_summary_stats_pvalue_gridF_scen3_n250$pvalue <- NA

min_mse_gridF_scen3_n250 <- min(df_summary_stats_gridF_scen3_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen3_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen3_n250[[id+1]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen3_n250 <- min_mse_gridF_scen3_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen3_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen3_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen3_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen3_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * covariate_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                            par = c(comb$mean_param1, comb$mean_param2), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen3_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen3_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen3_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen3_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen3_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen3_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen3_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen3_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen3_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen3_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen3_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen3_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen3_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen3_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen3_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen3_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen3_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen3_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 15 | nc = 10"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen3_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen3_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen3_n250)
dim(comb_subset)
length(which(comb_subset$pvalue >= 0.1))

save(df_summary_stats_pvalue_gridF_scen3_n250, file = "pvalue_plot_gridF_scen3_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen3_n250 <- df_summary_stats_gridT_scen3_n250
df_summary_stats_pvalue_gridT_scen3_n250$pvalue <- NA

min_mse_gridT_scen3_n250 <- min(df_summary_stats_gridT_scen3_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen3_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen3_n250[[id+1]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen3_n250 <- min_mse_gridT_scen3_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen3_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen3_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen3_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen3_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * covariate_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                            par = c(comb$mean_param1, comb$mean_param2), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen3_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen3_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen3_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen3_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen3_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen3_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen3_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen3_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen3_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen3_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen3_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen3_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen3_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen3_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen3_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen3_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen3_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen3_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 48 | nc = 19"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen3_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen3_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen3_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen3_n250, file = "pvalue_plot_gridT_scen3_n250.RData")





########################
### SCEN 3 - N = 500 ###
########################

true_params_scen3_n500
df_summary_stats_gridF_scen3_n500 <- df_summary_stats_scen3_n500[df_summary_stats_scen3_n500$grid == FALSE,]
df_summary_stats_gridT_scen3_n500 <- df_summary_stats_scen3_n500[df_summary_stats_scen3_n500$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen3_n500 <- df_summary_stats_gridF_scen3_n500
df_summary_stats_pvalue_gridF_scen3_n500$pvalue <- NA

min_mse_gridF_scen3_n500 <- min(df_summary_stats_gridF_scen3_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen3_n500)
count_significant <- 0
prop_explored <- 0   

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen3_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen3_n500 <- min_mse_gridF_scen3_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen3_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen3_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen3_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen3_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * covariate_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                            par = c(comb$mean_param1, comb$mean_param2), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen3_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen3_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen3_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen3_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen3_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen3_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen3_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen3_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen3_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen3_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen3_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen3_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen3_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen3_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen3_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen3_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen3_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen3_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 10 | nc = 11"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen3_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen3_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen3_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen3_n500, file = "pvalue_plot_gridF_scen3_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen3_n500 <- df_summary_stats_gridT_scen3_n500
df_summary_stats_pvalue_gridT_scen3_n500$pvalue <- NA

min_mse_gridT_scen3_n500 <- min(df_summary_stats_gridT_scen3_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen3_n500)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen3_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen3_n500 <- min_mse_gridT_scen3_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen3_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen3_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen3_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen3_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * covariate_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
                                            par = c(comb$mean_param1, comb$mean_param2), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen3_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen3_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen3_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen3_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen3_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen3_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen3_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen3_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen3_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen3_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen3_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen3_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen3_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen3_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen3_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen3_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen3_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen3_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 40 | nc = 17"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen3_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen3_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen3_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen3_n500, file = "pvalue_plot_gridT_scen3_n500.RData")





########################
### SCEN 4 - N = 100 ###
########################

true_params_scen4_n100
df_summary_stats_gridF_scen4_n100 <- df_summary_stats_scen4_n100[df_summary_stats_scen4_n100$grid == FALSE,]
df_summary_stats_gridT_scen4_n100 <- df_summary_stats_scen4_n100[df_summary_stats_scen4_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen4_n100 <- df_summary_stats_gridF_scen4_n100
df_summary_stats_pvalue_gridF_scen4_n100$pvalue <- NA

min_mse_gridF_scen4_n100 <- min(df_summary_stats_gridF_scen4_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen4_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen4_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen4_n100 <- min_mse_gridF_scen4_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen4_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen4_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen4_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen4_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * covariate1_values +
                            comb$mean_param3 * covariate2_values + comb$mean_param4 * covariate3_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])}, verbose = F,
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1 +
                        comb$mean_param3 * bestmod_ppp[[j]]$V2 + comb$mean_param4 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen4_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen4_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen4_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen4_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen4_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen4_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen4_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen4_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen4_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen4_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen4_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen4_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen4_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen4_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen4_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen4_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen4_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen4_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 26 | nc = 6"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen4_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen4_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen4_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen4_n100, file = "pvalue_plot_gridF_scen4_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen4_n100 <- df_summary_stats_gridT_scen4_n100
df_summary_stats_pvalue_gridT_scen4_n100$pvalue <- NA

min_mse_gridT_scen4_n100 <- min(df_summary_stats_gridT_scen4_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen4_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen4_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen4_n100 <- min_mse_gridT_scen4_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen4_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen4_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen4_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen4_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * covariate1_values +
                           comb$mean_param3 * covariate2_values + comb$mean_param4 * covariate3_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1 +
                        comb$mean_param3 * bestmod_ppp[[j]]$V2 + comb$mean_param4 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen4_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen4_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen4_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen4_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen4_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen4_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen4_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen4_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen4_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen4_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen4_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen4_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen4_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen4_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen4_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen4_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen4_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen4_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 4 | nc = 5"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen4_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen4_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen4_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen4_n100, file = "pvalue_plot_gridT_scen4_n100.RData")





########################
### SCEN 4 - N = 250 ###
########################

true_params_scen4_n250
df_summary_stats_gridF_scen4_n250 <- df_summary_stats_scen4_n250[df_summary_stats_scen4_n250$grid == FALSE,]
df_summary_stats_gridT_scen4_n250 <- df_summary_stats_scen4_n250[df_summary_stats_scen4_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen4_n250 <- df_summary_stats_gridF_scen4_n250
df_summary_stats_pvalue_gridF_scen4_n250$pvalue <- NA

min_mse_gridF_scen4_n250 <- min(df_summary_stats_gridF_scen4_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen4_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen4_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen4_n250 <- min_mse_gridF_scen4_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen4_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen4_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen4_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen4_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * covariate1_values +
                            comb$mean_param3 * covariate2_values + comb$mean_param4 * covariate3_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1 +
                        comb$mean_param3 * bestmod_ppp[[j]]$V2 + comb$mean_param4 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen4_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen4_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen4_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen4_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen4_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen4_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen4_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen4_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen4_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen4_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen4_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen4_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen4_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen4_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen4_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen4_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen4_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen4_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 41 | nc = 5"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen4_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen4_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen4_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen4_n250, file = "pvalue_plot_gridF_scen4_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen4_n250 <- df_summary_stats_gridT_scen4_n250
df_summary_stats_pvalue_gridT_scen4_n250$pvalue <- NA

min_mse_gridT_scen4_n250 <- min(df_summary_stats_gridT_scen4_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen4_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen4_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen4_n250 <- min_mse_gridT_scen4_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen4_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen4_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen4_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen4_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * covariate1_values +
                           comb$mean_param3 * covariate2_values + comb$mean_param4 * covariate3_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1 +
                        comb$mean_param3 * bestmod_ppp[[j]]$V2 + comb$mean_param4 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen4_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen4_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen4_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen4_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen4_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen4_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen4_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen4_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen4_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen4_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen4_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen4_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen4_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen4_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen4_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen4_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen4_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen4_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 25 | nc = 10"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen4_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen4_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen4_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen4_n250, file = "pvalue_plot_gridT_scen4_n250.RData")





########################
### SCEN 4 - N = 500 ###
########################

true_params_scen4_n500
df_summary_stats_gridF_scen4_n500 <- df_summary_stats_scen4_n500[df_summary_stats_scen4_n500$grid == FALSE,]
df_summary_stats_gridT_scen4_n500 <- df_summary_stats_scen4_n500[df_summary_stats_scen4_n500$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen4_n500 <- df_summary_stats_gridF_scen4_n500
df_summary_stats_pvalue_gridF_scen4_n500$pvalue <- NA

min_mse_gridF_scen4_n500 <- min(df_summary_stats_gridF_scen4_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen4_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen4_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen4_n500 <- min_mse_gridF_scen4_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen4_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen4_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen4_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen4_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * covariate1_values +
                            comb$mean_param3 * covariate2_values + comb$mean_param4 * covariate3_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1 +
                        comb$mean_param3 * bestmod_ppp[[j]]$V2 + comb$mean_param4 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen4_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen4_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen4_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen4_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen4_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen4_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen4_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen4_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen4_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen4_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen4_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen4_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen4_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen4_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen4_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen4_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen4_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen4_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 22 | nc = 7"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen4_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen4_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen4_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen4_n500, file = "pvalue_plot_gridF_scen4_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen4_n500 <- df_summary_stats_gridT_scen4_n500
df_summary_stats_pvalue_gridT_scen4_n500$pvalue <- NA

min_mse_gridT_scen4_n500 <- min(df_summary_stats_gridT_scen4_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen4_n500)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen4_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen4_n500 <- min_mse_gridT_scen4_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen4_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen4_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen4_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen4_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * covariate1_values +
                           comb$mean_param3 * covariate2_values + comb$mean_param4 * covariate3_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1] + a[3]*cov[2] + a[4]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$V1 +
                        comb$mean_param3 * bestmod_ppp[[j]]$V2 + comb$mean_param4 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen4_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen4_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen4_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen4_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen4_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen4_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen4_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen4_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen4_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen4_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen4_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen4_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen4_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen4_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen4_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen4_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen4_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen4_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 17 | nc = 11"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen4_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen4_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen4_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen4_n500, file = "pvalue_plot_gridT_scen4_n500.RData")





########################
### SCEN 5 - N = 100 ###
########################

true_params_scen5_n100
df_summary_stats_gridF_scen5_n100 <- df_summary_stats_scen5_n100[df_summary_stats_scen5_n100$grid == FALSE,]
df_summary_stats_gridT_scen5_n100 <- df_summary_stats_scen5_n100[df_summary_stats_scen5_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen5_n100 <- df_summary_stats_gridF_scen5_n100
df_summary_stats_pvalue_gridF_scen5_n100$pvalue <- NA

min_mse_gridF_scen5_n100 <- min(df_summary_stats_gridF_scen5_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen5_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen5_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen5_n100 <- min_mse_gridF_scen5_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen5_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen5_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen5_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen5_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x + comb$mean_param3 * false_ppp$df$y + comb$mean_param4 * false_ppp$df$t +
                            comb$mean_param5 * covariate1_values + comb$mean_param6 * covariate2_values + comb$mean_param7 * covariate3_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$x + comb$mean_param3 * bestmod_ppp[[j]]$y + comb$mean_param4 * bestmod_ppp[[j]]$t +
                        comb$mean_param5 * bestmod_ppp[[j]]$V1 + comb$mean_param6 * bestmod_ppp[[j]]$V2 + comb$mean_param7 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen5_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen5_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen5_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen5_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen5_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen5_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen5_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen5_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen5_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen5_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen5_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen5_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen5_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen5_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen5_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen5_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen5_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen5_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 3 | nc = 3"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen5_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen5_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen5_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen5_n100, file = "pvalue_plot_gridF_scen5_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen5_n100 <- df_summary_stats_gridT_scen5_n100
df_summary_stats_pvalue_gridT_scen5_n100$pvalue <- NA

min_mse_gridT_scen5_n100 <- min(df_summary_stats_gridT_scen5_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen5_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen5_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen5_n100 <- min_mse_gridT_scen5_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen5_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen5_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen5_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen5_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x + comb$mean_param3 * true_ppp$df$y + comb$mean_param4 * true_ppp$df$t + 
                           comb$mean_param5 * covariate1_values + comb$mean_param6 * covariate2_values + comb$mean_param7 * covariate3_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$x + comb$mean_param3 * bestmod_ppp[[j]]$y + comb$mean_param4 * bestmod_ppp[[j]]$t + 
                        comb$mean_param5 * bestmod_ppp[[j]]$V1 + comb$mean_param6 * bestmod_ppp[[j]]$V2 + comb$mean_param7 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen5_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen5_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen5_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen5_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen5_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen5_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen5_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen5_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen5_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen5_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen5_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen5_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen5_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen5_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen5_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen5_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen5_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen5_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 50 | nc = 13"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen5_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen5_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen5_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen5_n100, file = "pvalue_plot_gridT_scen5_n100.RData")





########################
### SCEN 5 - N = 250 ###
########################

true_params_scen5_n250
df_summary_stats_gridF_scen5_n250 <- df_summary_stats_scen5_n250[df_summary_stats_scen5_n250$grid == FALSE,]
df_summary_stats_gridT_scen5_n250 <- df_summary_stats_scen5_n250[df_summary_stats_scen5_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen5_n250 <- df_summary_stats_gridF_scen5_n250
df_summary_stats_pvalue_gridF_scen5_n250$pvalue <- NA

min_mse_gridF_scen5_n250 <- min(df_summary_stats_gridF_scen5_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen5_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen5_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen5_n250 <- min_mse_gridF_scen5_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen5_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen5_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen5_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen5_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x + comb$mean_param3 * false_ppp$df$y + comb$mean_param4 * false_ppp$df$t +
                            comb$mean_param5 * covariate1_values + comb$mean_param6 * covariate2_values + comb$mean_param7 * covariate3_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$x + comb$mean_param3 * bestmod_ppp[[j]]$y + comb$mean_param4 * bestmod_ppp[[j]]$t +
                        comb$mean_param5 * bestmod_ppp[[j]]$V1 + comb$mean_param6 * bestmod_ppp[[j]]$V2 + comb$mean_param7 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen5_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen5_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen5_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen5_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen5_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen5_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen5_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen5_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen5_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen5_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen5_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen5_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen5_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen5_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen5_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen5_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen5_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen5_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 4 | nc = 6"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen5_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen5_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen5_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen5_n250, file = "pvalue_plot_gridF_scen5_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen5_n250 <- df_summary_stats_gridT_scen5_n250
df_summary_stats_pvalue_gridT_scen5_n250$pvalue <- NA

min_mse_gridT_scen5_n250 <- min(df_summary_stats_gridT_scen5_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen5_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen5_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen5_n250 <- min_mse_gridT_scen5_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen5_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen5_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen5_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen5_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x + comb$mean_param3 * true_ppp$df$y + comb$mean_param4 * true_ppp$df$t + 
                           comb$mean_param5 * covariate1_values + comb$mean_param6 * covariate2_values + comb$mean_param7 * covariate3_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$x + comb$mean_param3 * bestmod_ppp[[j]]$y + comb$mean_param4 * bestmod_ppp[[j]]$t + 
                        comb$mean_param5 * bestmod_ppp[[j]]$V1 + comb$mean_param6 * bestmod_ppp[[j]]$V2 + comb$mean_param7 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen5_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen5_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen5_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen5_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen5_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen5_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen5_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen5_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen5_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen5_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen5_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen5_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen5_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen5_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen5_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen5_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen5_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen5_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 17 | nc = 13"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen5_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen5_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen5_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen5_n250, file = "pvalue_plot_gridT_scen5_n250.RData")





########################
### SCEN 5 - N = 500 ###
########################

true_params_scen5_n500
df_summary_stats_gridF_scen5_n500 <- df_summary_stats_scen5_n500[df_summary_stats_scen5_n500$grid == FALSE,]
df_summary_stats_gridT_scen5_n500 <- df_summary_stats_scen5_n500[df_summary_stats_scen5_n500$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen5_n500 <- df_summary_stats_gridF_scen5_n500
df_summary_stats_pvalue_gridF_scen5_n500$pvalue <- NA

min_mse_gridF_scen5_n500 <- min(df_summary_stats_gridF_scen5_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen5_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen5_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen5_n500 <- min_mse_gridF_scen5_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen5_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen5_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen5_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen5_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_false <- exp(comb$mean_param1 + comb$mean_param2 * false_ppp$df$x + comb$mean_param3 * false_ppp$df$y + comb$mean_param4 * false_ppp$df$t +
                            comb$mean_param5 * covariate1_values + comb$mean_param6 * covariate2_values + comb$mean_param7 * covariate3_values)
    Khat_false <- Khat_spatial.3D(false_ppp, lambda = pred_int_false, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$x + comb$mean_param3 * bestmod_ppp[[j]]$y + comb$mean_param4 * bestmod_ppp[[j]]$t +
                        comb$mean_param5 * bestmod_ppp[[j]]$V1 + comb$mean_param6 * bestmod_ppp[[j]]$V2 + comb$mean_param7 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen5_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen5_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen5_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen5_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen5_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen5_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen5_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen5_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen5_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen5_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen5_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen5_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen5_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen5_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen5_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen5_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen5_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen5_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 36 | nc = 14"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen5_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen5_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen5_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen5_n500, file = "pvalue_plot_gridF_scen5_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen5_n500 <- df_summary_stats_gridT_scen5_n500
df_summary_stats_pvalue_gridT_scen5_n500$pvalue <- NA

min_mse_gridT_scen5_n500 <- min(df_summary_stats_gridT_scen5_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen5_n500)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen5_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen5_n500 <- min_mse_gridT_scen5_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen5_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen5_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen5_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen5_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    covariate2_values <- interp3D_fast(df_point_obs, df_covs$cov2, p = 81, d = 3)
    covariate3_values <- interp3D_fast(df_point_obs, df_covs$cov3, p = 81, d = 3)
    pred_int_true <- exp(comb$mean_param1 + comb$mean_param2 * true_ppp$df$x + comb$mean_param3 * true_ppp$df$y + comb$mean_param4 * true_ppp$df$t + 
                           comb$mean_param5 * covariate1_values + comb$mean_param6 * covariate2_values + comb$mean_param7 * covariate3_values)
    Khat_true <- Khat_spatial.3D(true_ppp, lambda = pred_int_true, correction = "translate")
    set.seed(2)
    bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1] + a[6]*cov[2] + a[7]*cov[3])},
                                            par = c(comb$mean_param1, comb$mean_param2, comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 100, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      pred_int <- exp(comb$mean_param1 + comb$mean_param2 * bestmod_ppp[[j]]$x + comb$mean_param3 * bestmod_ppp[[j]]$y + comb$mean_param4 * bestmod_ppp[[j]]$t + 
                        comb$mean_param5 * bestmod_ppp[[j]]$V1 + comb$mean_param6 * bestmod_ppp[[j]]$V2 + comb$mean_param7 * bestmod_ppp[[j]]$V3)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen5_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen5_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen5_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen5_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen5_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen5_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen5_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen5_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen5_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen5_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen5_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen5_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen5_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen5_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen5_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen5_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen5_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen5_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 12 | nc = 15"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen5_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen5_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen5_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen5_n500, file = "pvalue_plot_gridT_scen5_n500.RData")





########################
### SCEN 6 - N = 100 ###
########################

true_params_scen6_n100
df_summary_stats_gridF_scen6_n100 <- df_summary_stats_scen6_n100[df_summary_stats_scen6_n100$grid == FALSE,]
df_summary_stats_gridT_scen6_n100 <- df_summary_stats_scen6_n100[df_summary_stats_scen6_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen6_n100 <- df_summary_stats_gridF_scen6_n100
df_summary_stats_pvalue_gridF_scen6_n100$pvalue <- NA

min_mse_gridF_scen6_n100 <- min(df_summary_stats_gridF_scen6_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen6_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen6_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen6_n100 <- min_mse_gridF_scen6_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen6_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen6_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen6_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen6_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec]
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param1, nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param2, nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param3, nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim]
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen6_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen6_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen6_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen6_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen6_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen6_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen6_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen6_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen6_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen6_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen6_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen6_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen6_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen6_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen6_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen6_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen6_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen6_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 43 | nc = 9"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen6_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen6_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen6_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen6_n100, file = "pvalue_plot_gridF_scen6_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen6_n100 <- df_summary_stats_gridT_scen6_n100
df_summary_stats_pvalue_gridT_scen6_n100$pvalue <- NA

min_mse_gridT_scen6_n100 <- min(df_summary_stats_gridT_scen6_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen6_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen6_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen6_n100 <- min_mse_gridT_scen6_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen6_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen6_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen6_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen6_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec]
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param1, nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param2, nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param3, nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim]
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen6_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen6_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen6_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen6_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen6_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen6_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen6_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen6_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen6_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen6_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen6_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen6_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen6_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen6_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen6_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen6_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen6_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen6_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 50 | nc = 14"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen6_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen6_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen6_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen6_n100, file = "pvalue_plot_gridT_scen6_n100.RData")





########################
### SCEN 6 - N = 250 ###
########################

true_params_scen6_n250
df_summary_stats_gridF_scen6_n250 <- df_summary_stats_scen6_n250[df_summary_stats_scen6_n250$grid == FALSE,]
df_summary_stats_gridT_scen6_n250 <- df_summary_stats_scen6_n250[df_summary_stats_scen6_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen6_n250 <- df_summary_stats_gridF_scen6_n250
df_summary_stats_pvalue_gridF_scen6_n250$pvalue <- NA

min_mse_gridF_scen6_n250 <- min(df_summary_stats_gridF_scen6_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen6_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen6_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen6_n250 <- min_mse_gridF_scen6_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen6_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen6_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen6_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen6_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec]
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param1, nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param2, nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param3, nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim]
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen6_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen6_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen6_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen6_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen6_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen6_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen6_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen6_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen6_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen6_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen6_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen6_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen6_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen6_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen6_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen6_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen6_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen6_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 29 | nc = 12"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen6_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen6_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen6_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen6_n250, file = "pvalue_plot_gridF_scen6_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen6_n250 <- df_summary_stats_gridT_scen6_n250
df_summary_stats_pvalue_gridT_scen6_n250$pvalue <- NA

min_mse_gridT_scen6_n250 <- min(df_summary_stats_gridT_scen6_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen6_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen6_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen6_n250 <- min_mse_gridT_scen6_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen6_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen6_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen6_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen6_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec]
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param1, nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param2, nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param3, nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim]
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen6_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen6_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen6_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen6_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen6_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen6_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen6_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen6_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen6_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen6_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen6_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen6_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen6_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen6_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen6_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen6_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen6_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen6_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 24 | nc = 15"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen6_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen6_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen6_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen6_n250, file = "pvalue_plot_gridT_scen6_n250.RData")





########################
### SCEN 6 - N = 500 ###
########################

true_params_scen6_n500
df_summary_stats_gridF_scen6_n500 <- df_summary_stats_scen6_n500[df_summary_stats_scen6_n500$grid == FALSE,]
df_summary_stats_gridT_scen6_n500 <- df_summary_stats_scen6_n500[df_summary_stats_scen6_n500$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen6_n500 <- df_summary_stats_gridF_scen6_n500
df_summary_stats_pvalue_gridF_scen6_n500$pvalue <- NA

min_mse_gridF_scen6_n500 <- min(df_summary_stats_gridF_scen6_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen6_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen6_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen6_n500 <- min_mse_gridF_scen6_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen6_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen6_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen6_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen6_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec]
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param1, nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param2, nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param3, nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim]
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen6_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen6_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen6_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen6_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen6_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen6_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen6_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen6_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen6_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen6_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen6_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen6_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen6_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen6_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen6_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen6_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen6_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen6_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 28 | nc = 15"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen6_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen6_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen6_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen6_n500, file = "pvalue_plot_gridF_scen6_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen6_n500 <- df_summary_stats_gridT_scen6_n500
df_summary_stats_pvalue_gridT_scen6_n500$pvalue <- NA

min_mse_gridT_scen6_n500 <- min(df_summary_stats_gridT_scen6_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen6_n500)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen6_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen6_n500 <- min_mse_gridT_scen6_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen6_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen6_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen6_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen6_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec]
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param1, nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param2, nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1])}, par = comb$mean_param3, nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim]
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen6_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen6_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen6_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen6_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen6_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen6_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen6_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen6_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen6_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen6_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen6_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen6_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen6_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen6_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen6_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen6_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen6_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen6_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 42 | nc = 27"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen6_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen6_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen6_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen6_n500, file = "pvalue_plot_gridT_scen6_n500.RData")





########################
### SCEN 7 - N = 100 ###
########################

true_params_scen7_n100
df_summary_stats_gridF_scen7_n100 <- df_summary_stats_scen7_n100[df_summary_stats_scen7_n100$grid == FALSE,]
df_summary_stats_gridT_scen7_n100 <- df_summary_stats_scen7_n100[df_summary_stats_scen7_n100$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen7_n100 <- df_summary_stats_gridF_scen7_n100
df_summary_stats_pvalue_gridF_scen7_n100$pvalue <- NA

min_mse_gridF_scen7_n100 <- min(df_summary_stats_gridF_scen7_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen7_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen7_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen7_n100 <- min_mse_gridF_scen7_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen7_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen7_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen7_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen7_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * false_ppp$df$x + comb$mean_param5 * false_ppp$df$y + comb$mean_param6 * false_ppp$df$t
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$df$x + comb$mean_param5 * bestmod_ppp[[j]]$df$y + comb$mean_param6 * bestmod_ppp[[j]]$df$t
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen7_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen7_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen7_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen7_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen7_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen7_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen7_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen7_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen7_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen7_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen7_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen7_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen7_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen7_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen7_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen7_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen7_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen7_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 49 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen7_n100

comb_subset <- df_summary_stats_pvalue_gridF_scen7_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen7_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen7_n100, file = "pvalue_plot_gridF_scen7_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen7_n100 <- df_summary_stats_gridT_scen7_n100
df_summary_stats_pvalue_gridT_scen7_n100$pvalue <- NA

min_mse_gridT_scen7_n100 <- min(df_summary_stats_gridT_scen7_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen7_n100)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen7_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen7_n100 <- min_mse_gridT_scen7_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen7_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen7_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen7_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen7_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * true_ppp$df$x + comb$mean_param5 * true_ppp$df$y + comb$mean_param6 * true_ppp$df$t
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$df$x + comb$mean_param5 * bestmod_ppp[[j]]$df$y + comb$mean_param6 * bestmod_ppp[[j]]$df$t
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen7_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen7_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen7_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen7_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen7_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen7_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen7_n100, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen7_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen7_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen7_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen7_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen7_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen7_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen7_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen7_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen7_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen7_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen7_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 46 | nc = 11"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen7_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen7_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen7_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen7_n100, file = "pvalue_plot_gridT_scen7_n100.RData")





########################
### SCEN 7 - N = 250 ###
########################

true_params_scen7_n250
df_summary_stats_gridF_scen7_n250 <- df_summary_stats_scen7_n250[df_summary_stats_scen7_n250$grid == FALSE,]
df_summary_stats_gridT_scen7_n250 <- df_summary_stats_scen7_n250[df_summary_stats_scen7_n250$grid == TRUE,]



### grid = F ###

df_summary_stats_pvalue_gridF_scen7_n250 <- df_summary_stats_gridF_scen7_n250
df_summary_stats_pvalue_gridF_scen7_n250$pvalue <- NA

min_mse_gridF_scen7_n250 <- min(df_summary_stats_gridF_scen7_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen7_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen7_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen7_n250 <- min_mse_gridF_scen7_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen7_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen7_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen7_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen7_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * false_ppp$df$x + comb$mean_param5 * false_ppp$df$y + comb$mean_param6 * false_ppp$df$t
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$df$x + comb$mean_param5 * bestmod_ppp[[j]]$df$y + comb$mean_param6 * bestmod_ppp[[j]]$df$t
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen7_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen7_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen7_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen7_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen7_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen7_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen7_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen7_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen7_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen7_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen7_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen7_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen7_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen7_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen7_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen7_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen7_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen7_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 42 | nc = 3"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen7_n250

comb_subset <- df_summary_stats_pvalue_gridF_scen7_n250 %>%
  filter(mse_combined <= threshold_mse_gridF_scen7_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen7_n250, file = "pvalue_plot_gridF_scen7_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen7_n250 <- df_summary_stats_gridT_scen7_n250
df_summary_stats_pvalue_gridT_scen7_n250$pvalue <- NA

min_mse_gridT_scen7_n250 <- min(df_summary_stats_gridT_scen7_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen7_n250)
count_significant <- 0
prop_explored <- 0  

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen7_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen7_n250 <- min_mse_gridT_scen7_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen7_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen7_n250 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen7_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen7_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * true_ppp$df$x + comb$mean_param5 * true_ppp$df$y + comb$mean_param6 * true_ppp$df$t
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$df$x + comb$mean_param5 * bestmod_ppp[[j]]$df$y + comb$mean_param6 * bestmod_ppp[[j]]$df$t
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen7_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen7_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen7_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen7_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen7_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen7_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen7_n250, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen7_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen7_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen7_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen7_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen7_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen7_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen7_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen7_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen7_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen7_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen7_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 12 | nc = 9"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen7_n250

comb_subset <- df_summary_stats_pvalue_gridT_scen7_n250 %>%
  filter(mse_combined <= threshold_mse_gridT_scen7_n250)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen7_n250, file = "pvalue_plot_gridT_scen7_n250.RData")





########################
### SCEN 7 - N = 500 ###
########################

true_params_scen7_n500
df_summary_stats_gridF_scen7_n500 <- df_summary_stats_scen7_n500[df_summary_stats_scen7_n500$grid == FALSE,]
df_summary_stats_gridT_scen7_n500 <- df_summary_stats_scen7_n500[df_summary_stats_scen7_n500$grid == TRUE,]

### grid = F ###

df_summary_stats_pvalue_gridF_scen7_n500 <- df_summary_stats_gridF_scen7_n500
df_summary_stats_pvalue_gridF_scen7_n500$pvalue <- NA

min_mse_gridF_scen7_n500 <- min(df_summary_stats_gridF_scen7_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen7_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen7_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen7_n500 <- min_mse_gridF_scen7_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen7_n500 %>%
    filter(mse_combined <= threshold_mse_gridF_scen7_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen7_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen7_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * false_ppp$df$x + comb$mean_param5 * false_ppp$df$y + comb$mean_param6 * false_ppp$df$t
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$df$x + comb$mean_param5 * bestmod_ppp[[j]]$df$y + comb$mean_param6 * bestmod_ppp[[j]]$df$t
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen7_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen7_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen7_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen7_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen7_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen7_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen7_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen7_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen7_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen7_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen7_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen7_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen7_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen7_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen7_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen7_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen7_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen7_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 49 | nc = 15"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen7_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen7_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen7_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen7_n500, file = "pvalue_plot_gridF_scen7_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen7_n500 <- df_summary_stats_gridT_scen7_n500
df_summary_stats_pvalue_gridT_scen7_n500$pvalue <- NA

min_mse_gridT_scen7_n500 <- min(df_summary_stats_gridT_scen7_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen7_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen7_n500[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen7_n500 <- min_mse_gridT_scen7_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen7_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen7_n500 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen7_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen7_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * true_ppp$df$x + comb$mean_param5 * true_ppp$df$y + comb$mean_param6 * true_ppp$df$t
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(2)
    g2_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    set.seed(3)
    g3_bestmod_ppp <- rstpp(lambda = function(x,y,t,a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t)}, par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6), nsim = 100)
    bestmod_ppp <- setNames(replicate(100, data.frame()), paste0("X", 1:100))
    for(si in 1:100) {
      bestmod_ppp[[si]] <- stpm(data.frame(x = c(g1_bestmod_ppp[[si]]$df$x, g2_bestmod_ppp[[si]]$df$x, g3_bestmod_ppp[[si]]$df$x),
                                           y = c(g1_bestmod_ppp[[si]]$df$y, g2_bestmod_ppp[[si]]$df$y, g3_bestmod_ppp[[si]]$df$y),
                                           t = c(g1_bestmod_ppp[[si]]$df$t, g2_bestmod_ppp[[si]]$df$t, g3_bestmod_ppp[[si]]$df$t),
                                           m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]]$df)), rep("b", nrow(g2_bestmod_ppp[[si]]$df)), rep("c", nrow(g3_bestmod_ppp[[si]]$df))))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)))
    for(j in 1:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$df$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$df$x + comb$mean_param5 * bestmod_ppp[[j]]$df$y + comb$mean_param6 * bestmod_ppp[[j]]$df$t
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]]$df[,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen7_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen7_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen7_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen7_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen7_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 20){
    message("Raggiunti 20 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen7_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen7_n500, 4), 
            ". Interrompo senza aver raggiunto 20 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 20 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen7_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen7_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen7_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen7_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen7_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen7_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen7_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen7_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen7_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen7_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen7_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 16 | nc = 16"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen7_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen7_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen7_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen7_n500, file = "pvalue_plot_gridT_scen7_n500.RData")





########################
### SCEN 8 - N = 100 ###
########################

true_params_scen8_n100
df_summary_stats_gridF_scen8_n100 <- df_summary_stats_scen8_n100[df_summary_stats_scen8_n100$grid == FALSE,]
df_summary_stats_gridT_scen8_n100 <- df_summary_stats_scen8_n100[df_summary_stats_scen8_n100$grid == TRUE,]

df_summary_stats_pvalue_gridF_scen8_n100[which(df_summary_stats_pvalue_gridF_scen8_n100$mult == 10 & df_summary_stats_pvalue_gridF_scen8_n100$ncube == 5),]
df_summary_stats_pvalue_gridF_scen8_n100[which(df_summary_stats_pvalue_gridF_scen8_n100$mult == 1 & df_summary_stats_pvalue_gridF_scen8_n100$ncube == 1),]

df_summary_stats_pvalue_gridF_scen8_n100[which(df_summary_stats_pvalue_gridF_scen8_n100$contour_class == "white"),]

### grid = F ###

df_summary_stats_pvalue_gridF_scen8_n100 <- df_summary_stats_gridF_scen8_n100
df_summary_stats_pvalue_gridF_scen8_n100$pvalue <- NA

min_mse_gridF_scen8_n100 <- min(df_summary_stats_gridF_scen8_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen8_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen8_n100[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen8_n100 <- min_mse_gridF_scen8_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen8_n100 %>%
    filter(mse_combined <= threshold_mse_gridF_scen8_n100 & is.na(pvalue) & ncube <= 18)
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen8_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen8_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * false_ppp$df$x + comb$mean_param5 * false_ppp$df$y + comb$mean_param6 * false_ppp$df$t + comb$mean_param7 * covariate1_values
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(2)
    g2_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(3)
    g3_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- setNames(replicate(101, data.frame()), paste0("X", 1:101))
    for(si in 1:101) {
      bestmod_ppp[[si]] <- data.frame(x = c(g1_bestmod_ppp[[si]]$x, g2_bestmod_ppp[[si]]$x, g3_bestmod_ppp[[si]]$x),
                                      y = c(g1_bestmod_ppp[[si]]$y, g2_bestmod_ppp[[si]]$y, g3_bestmod_ppp[[si]]$y),
                                      t = c(g1_bestmod_ppp[[si]]$t, g2_bestmod_ppp[[si]]$t, g3_bestmod_ppp[[si]]$t),
                                      cov1 = c(g1_bestmod_ppp[[si]]$V1, g2_bestmod_ppp[[si]]$V1, g3_bestmod_ppp[[si]]$V1),
                                      m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]])), rep("b", nrow(g2_bestmod_ppp[[si]])), rep("c", nrow(g3_bestmod_ppp[[si]])))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)-1))
    for(j in 2:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$x + comb$mean_param5 * bestmod_ppp[[j]]$y + comb$mean_param6 * bestmod_ppp[[j]]$t + comb$mean_param7 * bestmod_ppp[[j]]$cov1
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
      Khat_sim[j-1, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen8_n100$pvalue[
      df_summary_stats_pvalue_gridF_scen8_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen8_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen8_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen8_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 200){
    message("Raggiunti 200 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen8_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen8_n100, 4), 
            ". Interrompo senza aver raggiunto 200 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 200 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen8_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen8_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen8_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen8_n100 <- ggplot(df_summary_stats_pvalue_gridF_scen8_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen8_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen8_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen8_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen8_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen8_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen8_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 3 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen8_n100


df_summary_stats_pvalue_gridF_scen8_n100 %>% filter(!is.na(pvalue)) %>% summarize(max(mse_combined))
comb_subset <- df_summary_stats_pvalue_gridF_scen8_n100 %>%
  filter(mse_combined <= threshold_mse_gridF_scen8_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen8_n100, file = "pvalue_plot_gridF_scen8_n100.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen8_n100 <- df_summary_stats_gridT_scen8_n100
df_summary_stats_pvalue_gridT_scen8_n100$pvalue <- NA

min_mse_gridT_scen8_n100 <- min(df_summary_stats_gridT_scen8_n100$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen8_n100)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen8_n100[[id]]


while(count_significant < 100){
  threshold_mse_gridT_scen8_n100 <- min_mse_gridT_scen8_n100 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen8_n100 %>%
    filter(mse_combined <= threshold_mse_gridT_scen8_n100 & is.na(pvalue))
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen8_n100, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen8_n100 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * true_ppp$df$x + comb$mean_param5 * true_ppp$df$y + comb$mean_param6 * true_ppp$df$t + comb$mean_param7 * covariate1_values
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(2)
    g2_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(3)
    g3_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- setNames(replicate(101, data.frame()), paste0("X", 1:101))
    for(si in 1:101) {
      bestmod_ppp[[si]] <- data.frame(x = c(g1_bestmod_ppp[[si]]$x, g2_bestmod_ppp[[si]]$x, g3_bestmod_ppp[[si]]$x),
                                      y = c(g1_bestmod_ppp[[si]]$y, g2_bestmod_ppp[[si]]$y, g3_bestmod_ppp[[si]]$y),
                                      t = c(g1_bestmod_ppp[[si]]$t, g2_bestmod_ppp[[si]]$t, g3_bestmod_ppp[[si]]$t),
                                      cov1 = c(g1_bestmod_ppp[[si]]$V1, g2_bestmod_ppp[[si]]$V1, g3_bestmod_ppp[[si]]$V1),
                                      m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]])), rep("b", nrow(g2_bestmod_ppp[[si]])), rep("c", nrow(g3_bestmod_ppp[[si]])))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)-1))
    for(j in 2:length(bestmod_ppp)){
      marks_vec_sim <- bestmod_ppp[[j]]$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$x + comb$mean_param5 * bestmod_ppp[[j]]$y + comb$mean_param6 * bestmod_ppp[[j]]$t + comb$mean_param7 * bestmod_ppp[[j]]$cov1
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j-1, ] <- Khat_sim_j$Khat
    }
    curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
    global_env <- global_envelope_test(curve_set, alternative = "two.sided")
    pvalue <- attr(global_env, "p")
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen8_n100$pvalue[
      df_summary_stats_pvalue_gridT_scen8_n100$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen8_n100$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen8_n100$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen8_n100, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 100){
    message("Raggiunti 100 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen8_n100, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 100){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen8_n100, 4), 
            ". Interrompo senza aver raggiunto 100 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 100 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen8_n100 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen8_n100$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen8_n100$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen8_n100 <- ggplot(df_summary_stats_pvalue_gridT_scen8_n100, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen8_n100, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen8_n100$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen8_n100$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen8_n100, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen8_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen8_n100, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 10 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen8_n100

comb_subset <- df_summary_stats_pvalue_gridT_scen8_n100 %>%
  filter(mse_combined <= threshold_mse_gridT_scen8_n100)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen8_n100, file = "pvalue_plot_gridT_scen8_n100.RData")





########################
### SCEN 8 - N = 250 ###
########################

true_params_scen8_n250
df_summary_stats_gridF_scen8_n250 <- df_summary_stats_scen8_n250[df_summary_stats_scen8_n250$grid == FALSE,]
df_summary_stats_gridT_scen8_n250 <- df_summary_stats_scen8_n250[df_summary_stats_scen8_n250$grid == TRUE,]

### grid = F ###

df_summary_stats_pvalue_gridF_scen8_n250 <- df_summary_stats_gridF_scen8_n250
df_summary_stats_pvalue_gridF_scen8_n250$pvalue <- NA

min_mse_gridF_scen8_n250 <- min(df_summary_stats_gridF_scen8_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen8_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen8_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridF_scen8_n250 <- min_mse_gridF_scen8_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen8_n250 %>%
    filter(mse_combined <= threshold_mse_gridF_scen8_n250 & is.na(pvalue) & ncube <= 18)
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen8_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen8_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- false_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * false_ppp$df$x + comb$mean_param5 * false_ppp$df$y + comb$mean_param6 * false_ppp$df$t + comb$mean_param7 * covariate1_values
    pred_int_false <- as.numeric(exp(linear_predictor))
    Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(2)
    g2_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(3)
    g3_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- setNames(replicate(101, data.frame()), paste0("X", 1:101))
    for(si in 1:101) {
      bestmod_ppp[[si]] <- data.frame(x = c(g1_bestmod_ppp[[si]]$x, g2_bestmod_ppp[[si]]$x, g3_bestmod_ppp[[si]]$x),
                                      y = c(g1_bestmod_ppp[[si]]$y, g2_bestmod_ppp[[si]]$y, g3_bestmod_ppp[[si]]$y),
                                      t = c(g1_bestmod_ppp[[si]]$t, g2_bestmod_ppp[[si]]$t, g3_bestmod_ppp[[si]]$t),
                                      cov1 = c(g1_bestmod_ppp[[si]]$V1, g2_bestmod_ppp[[si]]$V1, g3_bestmod_ppp[[si]]$V1),
                                      m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]])), rep("b", nrow(g2_bestmod_ppp[[si]])), rep("c", nrow(g3_bestmod_ppp[[si]])))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)-1))
    for (j in 2:length(bestmod_ppp)) {
      tryCatch({marks_vec_sim <- bestmod_ppp[[j]]$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] +comb$mean_param4 * bestmod_ppp[[j]]$x +comb$mean_param5 * bestmod_ppp[[j]]$y +comb$mean_param6 * bestmod_ppp[[j]]$t +comb$mean_param7 * bestmod_ppp[[j]]$cov1
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][, c(1,2,3)]),lambda = pred_int,correction = "translate",dist = Khat_false$dist)
      Khat_sim[j-1, ] <- Khat_sim_j$Khat
      },
      error = function(e) {message(paste("Errore in iterazione", j, ":", e$message)) 
        Khat_sim[j-1, ] <- NA
      })
    }
    pvalue <- tryCatch({
      curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
      global_env <- global_envelope_test(curve_set, alternative = "two.sided")
      attr(global_env, "p")
    },
    error = function(e) {
      message("Errore in global_envelope_test: assegno pvalue = 0.00989")
      0.00989
    })
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen8_n250$pvalue[
      df_summary_stats_pvalue_gridF_scen8_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen8_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen8_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen8_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 200){
    message("Raggiunti 200 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen8_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen8_n250, 4), 
            ". Interrompo senza aver raggiunto 200 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 200 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen8_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen8_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen8_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

combs_to_na <- data.frame(
  mult = c(rep(1,7), rep(2,6), rep(3,5), rep(4,4), rep(5,3), 6,6,7),
  ncube = c(11:17,12:17, 13:17, 14:17, 15:17,16,17,17)
)
df_summary_stats_pvalue_gridF_scen8_n250_mod <- set_contour_class_na(df_summary_stats_pvalue_gridF_scen8_n250_mod, combs_to_na)

rows_na <- df_summary_stats_pvalue_gridF_scen8_n250_mod[!is.na(df_summary_stats_pvalue_gridF_scen8_n250_mod$pvalue), ]

# Trova il valore massimo di mse_combined tra queste righe
max_mse <- max(rows_na$mse_combined, na.rm = TRUE)

max_mse

p_gridF_scen8_n250 <- ggplot(df_summary_stats_pvalue_gridF_scen8_n250_mod, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen8_n250_mod, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen8_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen8_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen8_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen8_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(1.013, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 3 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen8_n250


df_summary_stats_pvalue_gridF_scen8_n250_mod %>% filter(!is.na(contour_class)) %>% summarize(max(mse_combined))

comb_subset <- df_summary_stats_pvalue_gridF_scen8_n250_mod %>%
  filter(!is.na(contour_class))
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen8_n250, file = "pvalue_plot_gridF_scen8_n250.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen8_n250 <- df_summary_stats_gridT_scen8_n250
df_summary_stats_pvalue_gridT_scen8_n250$pvalue <- NA

min_mse_gridT_scen8_n250 <- min(df_summary_stats_gridT_scen8_n250$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen8_n250)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen8_n250[[id]]


while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen8_n250 <- min_mse_gridT_scen8_n250 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen8_n250 %>%
    filter(mse_combined <= threshold_mse_gridT_scen8_n250 & is.na(pvalue) & ncube <= 20)
  
  if(nrow(comb_subset) == 0){
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen8_n250, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen8_n250 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), ")")
    next
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * true_ppp$df$x + comb$mean_param5 * true_ppp$df$y + comb$mean_param6 * true_ppp$df$t + comb$mean_param7 * covariate1_values
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(2)
    g2_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(3)
    g3_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- setNames(replicate(101, data.frame()), paste0("X", 1:101))
    for(si in 1:101) {
      bestmod_ppp[[si]] <- data.frame(x = c(g1_bestmod_ppp[[si]]$x, g2_bestmod_ppp[[si]]$x, g3_bestmod_ppp[[si]]$x),
                                      y = c(g1_bestmod_ppp[[si]]$y, g2_bestmod_ppp[[si]]$y, g3_bestmod_ppp[[si]]$y),
                                      t = c(g1_bestmod_ppp[[si]]$t, g2_bestmod_ppp[[si]]$t, g3_bestmod_ppp[[si]]$t),
                                      cov1 = c(g1_bestmod_ppp[[si]]$V1, g2_bestmod_ppp[[si]]$V1, g3_bestmod_ppp[[si]]$V1),
                                      m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]])), rep("b", nrow(g2_bestmod_ppp[[si]])), rep("c", nrow(g3_bestmod_ppp[[si]])))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)-1))
    for(j in 2:length(bestmod_ppp)){
      tryCatch({ marks_vec_sim <- bestmod_ppp[[j]]$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$x + comb$mean_param5 * bestmod_ppp[[j]]$y + comb$mean_param6 * bestmod_ppp[[j]]$t + comb$mean_param7 * bestmod_ppp[[j]]$cov1
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j-1, ] <- Khat_sim_j$Khat
      },
      error = function(e) {message(paste("Errore in iterazione", j, ":", e$message)) 
        Khat_sim[j-1, ] <- NA
      })
    }
    pvalue <- tryCatch({
      curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
      global_env <- global_envelope_test(curve_set, alternative = "two.sided")
      attr(global_env, "p")
    },
    error = function(e) {
      message("Errore in global_envelope_test: assegno pvalue = 0.00989")
      0.00989
    })
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen8_n250$pvalue[
      df_summary_stats_pvalue_gridT_scen8_n250$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen8_n250$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen8_n250$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen8_n250, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 200){
    message("Raggiunti 100 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen8_n250, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen8_n250, 4), 
            ". Interrompo senza aver raggiunto 200 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 200 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen8_n250 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen8_n250$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen8_n250$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))


p_gridT_scen8_n250 <- ggplot(df_summary_stats_pvalue_gridT_scen8_n250, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen8_n250, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen8_n250$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen8_n250$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen8_n250, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen8_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen8_n250, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 9 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen8_n250

df_summary_stats_pvalue_gridT_scen8_n250 %>% filter(!is.na(contour_class)) %>% summarize(max(mse_combined))

comb_subset <- df_summary_stats_pvalue_gridT_scen8_n250 %>%
  filter(!is.na(contour_class))
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen8_n250, file = "pvalue_plot_gridT_scen8_n250.RData")





########################
### SCEN 8 - N = 500 ###
########################

true_params_scen8_n500
df_summary_stats_gridF_scen8_n500 <- df_summary_stats_scen8_n500[df_summary_stats_scen8_n500$grid == FALSE,]
df_summary_stats_gridT_scen8_n500 <- df_summary_stats_scen8_n500[df_summary_stats_scen8_n500$grid == TRUE,]

### grid = F ###

df_summary_stats_pvalue_gridF_scen8_n500 <- df_summary_stats_gridF_scen8_n500
df_summary_stats_pvalue_gridF_scen8_n500$pvalue <- NA

min_mse_gridF_scen8_n500 <- min(df_summary_stats_gridF_scen8_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridF_scen8_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
false_ppp <- inh_scen8_n500[[id]]

no_comb_counter <- 0

while(count_significant < 20 || prop_explored < 0.23){
  threshold_mse_gridF_scen8_n500 <- min_mse_gridF_scen8_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridF_scen8_n500 %>%
    filter(ncube == 8 & mult == 41)
  
  if(nrow(comb_subset) == 0){
    no_comb_counter <- no_comb_counter + 1
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridF_scen8_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridF_scen8_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), 
            ") | Tentativi consecutivi senza combinazioni = ", no_comb_counter)
    
    if(no_comb_counter >= 50){
      message("⚠️  Nessuna combinazione trovata per 50 aumenti consecutivi del threshold. Interrompo il ciclo.")
      break
    }
    next
  } else {
    no_comb_counter <- 0  # reset contatore se troviamo combinazioni
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind,
                     .packages = c("spatstat", "stpp", "stopp", "mgcv", "GET"),
                     .export   = c("df_covs", "interp3D_fast", "rstpp_cov_generalised_v2",
                                   "Khat_spatial.3D")) %dopar% {
                                     
                                     comb <- comb_subset[i, ]
                                     df_point_obs <- data.frame(x = false_ppp$df$x, y = false_ppp$df$y, t = false_ppp$df$t)
                                     covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
                                     param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
                                     marks_vec <- false_ppp$df$m1
                                     linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * false_ppp$df$x + comb$mean_param5 * false_ppp$df$y + comb$mean_param6 * false_ppp$df$t + comb$mean_param7 * covariate1_values
                                     pred_int_false <- as.numeric(exp(linear_predictor))
                                     Khat_false <- Khat_spatial.3D(stp(false_ppp$df[,c(1,2,3)]), lambda = pred_int_false, correction = "translate")
                                     set.seed(1)
                                     g1_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                                                par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
                                     set.seed(2)
                                     g2_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                                                par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
                                     set.seed(3)
                                     g3_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                                                                par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
                                     bestmod_ppp <- setNames(replicate(101, data.frame()), paste0("X", 1:101))
                                     for(si in 1:101) {
                                       bestmod_ppp[[si]] <- data.frame(x = c(g1_bestmod_ppp[[si]]$x, g2_bestmod_ppp[[si]]$x, g3_bestmod_ppp[[si]]$x),
                                                                       y = c(g1_bestmod_ppp[[si]]$y, g2_bestmod_ppp[[si]]$y, g3_bestmod_ppp[[si]]$y),
                                                                       t = c(g1_bestmod_ppp[[si]]$t, g2_bestmod_ppp[[si]]$t, g3_bestmod_ppp[[si]]$t),
                                                                       cov1 = c(g1_bestmod_ppp[[si]]$V1, g2_bestmod_ppp[[si]]$V1, g3_bestmod_ppp[[si]]$V1),
                                                                       m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]])), rep("b", nrow(g2_bestmod_ppp[[si]])), rep("c", nrow(g3_bestmod_ppp[[si]])))))
                                     }
                                     #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
                                     Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_false$dist), nrow = length(bestmod_ppp)-1))
                                     for(j in 2:length(bestmod_ppp)){
                                       tryCatch({ marks_vec_sim <- bestmod_ppp[[j]]$m1
                                       linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$x + comb$mean_param5 * bestmod_ppp[[j]]$y + comb$mean_param6 * bestmod_ppp[[j]]$t + comb$mean_param7 * bestmod_ppp[[j]]$cov1
                                       pred_int <- exp(linear_predictor_sim)
                                       Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_false$dist)
                                       Khat_sim[j-1, ] <- Khat_sim_j$Khat
                                       },
                                       error = function(e) {message(paste("Errore in iterazione", j, ":", e$message)) 
                                         Khat_sim[j-1, ] <- NA
                                       })
                                     }
                                     pvalue <- tryCatch({
                                       curve_set <- create_curve_set(list(r = Khat_false$dist, obs = Khat_false$Khat, sim_m = t(Khat_sim)))
                                       global_env <- global_envelope_test(curve_set, alternative = "two.sided")
                                       attr(global_env, "p")
                                     },
                                     error = function(e) {
                                       message("Errore in global_envelope_test: assegno pvalue = 0.00989")
                                       0.00989
                                     })
                                     c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
                                   }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridF_scen8_n500$pvalue[
      df_summary_stats_pvalue_gridF_scen8_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridF_scen8_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridF_scen8_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridF_scen8_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 203){
    message("Raggiunti 203 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridF_scen8_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.23){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridF_scen8_n500, 4), 
            ". Interrompo senza aver raggiunto 203 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 203 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridF_scen8_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridF_scen8_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridF_scen8_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridF_scen8_n500 <- ggplot(df_summary_stats_pvalue_gridF_scen8_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridF_scen8_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen8_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen8_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_false_scen8_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridF_scen8_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridF_scen8_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 3 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridF_scen8_n500

comb_subset <- df_summary_stats_pvalue_gridF_scen8_n500 %>%
  filter(mse_combined <= threshold_mse_gridF_scen8_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridF_scen8_n500, file = "pvalue_plot_gridF_scen8_n500.RData")



### grid = T ###

df_summary_stats_pvalue_gridT_scen8_n500 <- df_summary_stats_gridT_scen8_n500
df_summary_stats_pvalue_gridT_scen8_n500$pvalue <- NA

min_mse_gridT_scen8_n500 <- min(df_summary_stats_gridT_scen8_n500$mse_combined)
threshold_factor <- 1.10  
total_n <- nrow(df_summary_stats_pvalue_gridT_scen8_n500)
count_significant <- 0
prop_explored <- 0 

set.seed(2)
id <- sample(1:100, size = 1)
true_ppp <- inh_scen8_n500[[id]]

no_comb_counter <- 0

while(count_significant < 20 || prop_explored < 0.2){
  threshold_mse_gridT_scen8_n500 <- min_mse_gridT_scen8_n500 * threshold_factor
  comb_subset <- df_summary_stats_pvalue_gridT_scen8_n500 %>%
    filter(mse_combined <= threshold_mse_gridT_scen8_n500 & is.na(pvalue) & ncube <= 30)
  
  if(nrow(comb_subset) == 0){
    no_comb_counter <- no_comb_counter + 1
    threshold_factor <- threshold_factor + 0.05
    message("Nessuna combinazione entro la soglia ", round(threshold_mse_gridT_scen8_n500, 4), 
            ". Incremento soglia del 5%: nuova soglia = ", round(min_mse_gridT_scen8_n500 * threshold_factor, 4),
            " (fattore = ", round(threshold_factor, 3), 
            ") | Tentativi consecutivi senza combinazioni = ", no_comb_counter)
    
    if(no_comb_counter >= 50){
      message("⚠️  Nessuna combinazione trovata per 50 aumenti consecutivi del threshold. Interrompo il ciclo.")
      break
    }
    next
  } else {
    no_comb_counter <- 0  # reset contatore se troviamo combinazioni
  }
  
  pb <- progress_bar$new(
    format = " [:bar] :percent | Iterazione :current/:total | Tempo trascorso: :elapsed | ETA: :eta", 
    total = nrow(comb_subset), clear = FALSE, width = 80)
  
  results <- foreach(i = 1:nrow(comb_subset), .combine = rbind, .packages = c("spatstat", "stpp", "stopp", "mgcv")) %dopar% {
    
    comb <- comb_subset[i, ]
    df_point_obs <- data.frame(x = true_ppp$df$x, y = true_ppp$df$y, t = true_ppp$df$t)
    covariate1_values <- interp3D_fast(df_point_obs, df_covs$cov1, p = 81, d = 3)
    param_vec <- c(a = comb$mean_param1, b = comb$mean_param2, c = comb$mean_param3)
    marks_vec <- true_ppp$df$m1
    linear_predictor <- param_vec[marks_vec] + comb$mean_param4 * true_ppp$df$x + comb$mean_param5 * true_ppp$df$y + comb$mean_param6 * true_ppp$df$t + comb$mean_param7 * covariate1_values
    pred_int_true <- as.numeric(exp(linear_predictor))
    Khat_true <- Khat_spatial.3D(stp(true_ppp$df[,c(1,2,3)]), lambda = pred_int_true, correction = "translate")
    set.seed(1)
    g1_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param1, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(2)
    g2_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param2, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    set.seed(3)
    g3_bestmod_ppp <- rstpp_cov_generalised_v2(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])},
                                               par = c(comb$mean_param3, comb$mean_param4, comb$mean_param5, comb$mean_param6, comb$mean_param7), nsim = 101, covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))
    bestmod_ppp <- setNames(replicate(101, data.frame()), paste0("X", 1:101))
    for(si in 1:101) {
      bestmod_ppp[[si]] <- data.frame(x = c(g1_bestmod_ppp[[si]]$x, g2_bestmod_ppp[[si]]$x, g3_bestmod_ppp[[si]]$x),
                                      y = c(g1_bestmod_ppp[[si]]$y, g2_bestmod_ppp[[si]]$y, g3_bestmod_ppp[[si]]$y),
                                      t = c(g1_bestmod_ppp[[si]]$t, g2_bestmod_ppp[[si]]$t, g3_bestmod_ppp[[si]]$t),
                                      cov1 = c(g1_bestmod_ppp[[si]]$V1, g2_bestmod_ppp[[si]]$V1, g3_bestmod_ppp[[si]]$V1),
                                      m1 = as.factor(c(rep("a", nrow(g1_bestmod_ppp[[si]])), rep("b", nrow(g2_bestmod_ppp[[si]])), rep("c", nrow(g3_bestmod_ppp[[si]])))))
    }
    #bestmod_ppp <- Filter(function(df) all(complete.cases(df)), bestmod_ppp)
    Khat_sim <- as.data.frame(matrix(NA, ncol = length(Khat_true$dist), nrow = length(bestmod_ppp)-1))
    for(j in 2:length(bestmod_ppp)){
      tryCatch({ marks_vec_sim <- bestmod_ppp[[j]]$m1
      linear_predictor_sim <- param_vec[marks_vec_sim] + comb$mean_param4 * bestmod_ppp[[j]]$x + comb$mean_param5 * bestmod_ppp[[j]]$y + comb$mean_param6 * bestmod_ppp[[j]]$t + comb$mean_param7 * bestmod_ppp[[j]]$cov1
      pred_int <- exp(linear_predictor_sim)
      Khat_sim_j <- Khat_spatial.3D(stp(bestmod_ppp[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_true$dist)
      Khat_sim[j-1, ] <- Khat_sim_j$Khat
      },
      error = function(e) {message(paste("Errore in iterazione", j, ":", e$message)) 
        Khat_sim[j-1, ] <- NA
      })
    }
    pvalue <- tryCatch({
      curve_set <- create_curve_set(list(r = Khat_true$dist, obs = Khat_true$Khat, sim_m = t(Khat_sim)))
      global_env <- global_envelope_test(curve_set, alternative = "two.sided")
      attr(global_env, "p")
    },
    error = function(e) {
      message("Errore in global_envelope_test: assegno pvalue = 0.00989")
      0.00989
    })
    c(mult = comb$mult, ncube = comb$ncube, pvalue = pvalue)
  }
  if (is.null(dim(results))) {results <- as.data.frame(t(results))}
  for(r in 1:nrow(results)){
    df_summary_stats_pvalue_gridT_scen8_n500$pvalue[
      df_summary_stats_pvalue_gridT_scen8_n500$mult == results[r, "mult"] &
        df_summary_stats_pvalue_gridT_scen8_n500$ncube == results[r, "ncube"]
    ] <- results[r, "pvalue"]
  }
  new_significant <- sum(!is.na(results[, "pvalue"]) & results[, "pvalue"] <= 0.01)
  count_significant <- count_significant + new_significant
  prop_explored <- mean(!is.na(df_summary_stats_pvalue_gridT_scen8_n500$pvalue))
  
  message("Esplorato finora: ", round(100 * prop_explored, 1), "% del dataset")
  
  message("Trovati ", new_significant, " nuovi p-value ≤ 0.01 (totale finora: ", count_significant, 
          ") | Soglia MSE = ", round(threshold_mse_gridT_scen8_n500, 4), " | Fattore = ", round(threshold_factor, 3))
  
  if(count_significant >= 200){
    message("Raggiunti 200 p-value ≤ 0.01 entro la soglia ", round(threshold_mse_gridT_scen8_n500, 4), 
            " (fattore = ", round(threshold_factor, 3), ")")
    break
  }
  if(nrow(comb_subset) == 0 && count_significant < 20 || prop_explored < 0.2){
    message("Nessuna combinazione rimasta entro la soglia ", round(threshold_mse_gridT_scen8_n500, 4), 
            ". Interrompo senza aver raggiunto 200 p-value ≤ 0.01.")
    break
  }
  threshold_factor <- threshold_factor + 0.05
  message("Ancora meno di 200 p-value ≤ 0.01. Nuova soglia = ", round(min_mse_gridT_scen8_n500 * threshold_factor, 4), 
          " (fattore = ", round(threshold_factor, 3), ")")
}


df_summary_stats_pvalue_gridT_scen8_n500$contour_class <- cut(
  df_summary_stats_pvalue_gridT_scen8_n500$pvalue,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("red", "orange", "yellow", "white"))

p_gridT_scen8_n500 <- ggplot(df_summary_stats_pvalue_gridT_scen8_n500, aes(x = mult, y = ncube, fill = mse_combined)) +
  geom_tile() +
  geom_tile(data = subset(df_summary_stats_pvalue_gridT_scen8_n500, !is.na(contour_class)),
            aes(color = contour_class), fill = NA, linewidth = 1, width = 1, height = 1) +
  scale_fill_gradient(low = "blue", high = "brown4", na.value = "grey90",
                      limits = c(0, max(df_summary_stats_scen8_n500$mse_combined)),
                      breaks = seq(0, max(df_summary_stats_scen8_n500$mse_combined), length.out = 8),
                      labels = function(x) round(x),
                      guide = guide_colorbar(title = "MSE", title.theme = element_text(size = 16), label.theme = element_text(size = 14), order = 2)) +
  scale_color_manual(values = c("white" = "white", "yellow" = "yellow", "orange" = "orange", "red" = "red"),
                     name   = "p-value", breaks = c("white", "yellow", "orange", "red"), labels = c("p \u2265 0.10", "0.05 \u2264 p < 0.10", "0.01 \u2264 p < 0.05", "p < 0.01"),
                     guide  = guide_legend(order = 1, title.theme = element_text(size = 18), label.theme = element_text(size = 16), override.aes = list(fill   = c("white", "yellow", "orange", "red"), colour = "black", size   = 5), keywidth  = grid::unit(0.9, "cm"), keyheight = grid::unit(0.9, "cm"))) +
  geom_point(data = min_point_true_scen8_n500, aes(x = mult, y = ncube),
             shape = 21, fill = "green", color = "black", size = 5, stroke = 1, inherit.aes = FALSE) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c])) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(5, "pt"),
        legend.title = element_text(size = 16),    
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.position = "right",
        legend.box = "vertical") +
  annotate("text", x = 48, y = 48, label = bquote(MSE[min] == .(round(min_mse_gridT_scen8_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 45, label = bquote(MSE[max] == .(round(threshold_mse_gridT_scen8_n500, 3))), hjust = 1, size = 7, color = "white") +
  annotate("text", x = 48, y = 42, label = paste("q = 8 | nc = 1"), hjust = 1, size = 7, color = "green", fontface = "bold") +
  scale_x_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5), sec.axis = dup_axis(name = NULL, labels = NULL)) + 
  coord_fixed()

p_gridT_scen8_n500

comb_subset <- df_summary_stats_pvalue_gridT_scen8_n500 %>%
  filter(mse_combined <= threshold_mse_gridT_scen8_n500)
dim(comb_subset)[1]/2500
length(which(comb_subset$pvalue >= 0.1))/2500

save(df_summary_stats_pvalue_gridT_scen8_n500, file = "pvalue_plot_gridT_scen8_n500.RData")
