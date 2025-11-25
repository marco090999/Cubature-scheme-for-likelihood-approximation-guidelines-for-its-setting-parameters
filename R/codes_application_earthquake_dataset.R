

###################################################
#### REAL DATA APPLICATION: EARTHQUAKE DATASET ####
###################################################

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(plotly)
library(scales)

ppp_appl_terr <- terremoti_stp[which(pb > 0.99),]

# --- Prep: rename columns for clarity ---
quakes <- ppp_appl_terr %>%
  rename(lon = long, lat = lat, depth = z, mag = magn1) %>%
  filter(is.finite(lon), is.finite(lat), is.finite(depth), is.finite(mag))

faults <- faglie %>%
  rename(lon = x, lat = y) %>%
  filter(is.finite(lon), is.finite(lat))

# --- Convert to sf (WGS84) ---
quakes_sf <- st_as_sf(quakes, coords = c("lon","lat"), crs = 4326, remove = FALSE)
faults_sf <- st_as_sf(faults, coords = c("lon","lat"), crs = 4326, remove = FALSE)

# --- Greece polygon (for basemap) ---
gr <- ne_countries(country = "Greece", scale = "medium", returnclass = "sf")

# --- Map extent from data (+ small padding) ---
bbox_combined <- st_bbox(st_union(st_geometry(quakes_sf), st_geometry(faults_sf)))
pad_x <- 0.5; pad_y <- 0.5
xlim <- c(as.numeric(bbox_combined["xmin"]) - pad_x,
          as.numeric(bbox_combined["xmax"]) + pad_x)
ylim <- c(as.numeric(bbox_combined["ymin"]) - pad_y,
          as.numeric(bbox_combined["ymax"]) + pad_y)

# ==========================
# 1) 2D MAP: earthquakes + faults + Greece basemap
#    - x = lon, y = lat
#    - earthquakes colored by magnitude (continuous)
#    - faults in blue
# ==========================

mag_range <- range(quakes$mag, na.rm = TRUE)

p_map <- ggplot() +
  geom_sf(data = gr, fill = "grey95", color = "grey70", linewidth = 0.4) +
  geom_point(data = faults, aes(x = lon, y = lat),
             color = "#1f77b4", size = 0.1, alpha = 0.8) +
  geom_point(data = quakes, aes(x = lon, y = lat, color = mag),
             size = 1.6, alpha = 0.9) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  scale_color_gradient(
    low = "orange", high = "darkred",
    name = "Magnitude", limits = mag_range, breaks = seq(3,6.5, by =0.5),
    guide = guide_colorbar(
      barheight = unit(2.5, "npc"),  # ~ altezza pannello
      barwidth  = unit(16, "pt"),
      ticks.colour = "black"
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 18) +
  theme(
    # legenda/colorbar alta e aderente
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    
    # bordo nero spesso su tutti i lati (linee top/right senza etichette)
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    
    # tick spessi su assi principali (bottom/left)
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(5, "pt"),
    
    # niente etichette su assi top/right (frame resta visibile)
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank()
  )

p_map


# ==========================
# 2) 3D SCATTER: x = lon, y = lat, z = depth
#    - color by magnitude treated as categorical
#    - breaks = c(3, 3.2, 4, 7)
#    - depth axis reversed so 'down' is deeper
# ==========================

quakes <- ppp_appl_terr %>%
  rename(lon = long, lat = lat, depth = z, mag = magn1) %>%
  filter(is.finite(lon), is.finite(lat), is.finite(depth), is.finite(mag))

quakes_cat <- quakes %>%
  mutate(mag_cat = cut(mag,
                       breaks = c(3, 3.2, 4, 7),
                       include.lowest = TRUE, right = FALSE,
                       labels = c("[3, 3.2)", "[3.2, 4)", "[4, 7]")))
lvl  <- levels(quakes_cat$mag_cat)
cols <- c("orange", "red", "darkred"); names(cols) <- lvl
idx  <- as.integer(quakes_cat$mag_cat)
K    <- nlevels(quakes_cat$mag_cat)

op <- par(no.readonly = TRUE)
par(mar = c(3.5, 3.8, 2.2, 2.0),   # bottom, left, top, right (più piccoli => area plot più grande)
    mgp = c(2.2, 0.7, 0))          # distanze titolo/etichette/asse

## 3) PLOT 3D (ingrandisci etichette assi e tick)
scatter3D(
  x = quakes_cat$lon,
  y = quakes_cat$lat,
  z = quakes_cat$depth,       # usa -depth se vuoi profondità verso il basso
  theta = -45, phi = 40,
  pch = 20, cex = 0.8,        # punti un po' più grandi (se vuoi)
  ticktype = "detailed",
  colvar = idx, col = unname(cols), clim = c(1, K),
  colkey = FALSE,
  xlab = "x", ylab = "y", zlab = "z",
  bty  = "f", lwd = 1.8,      # box più spesso
  cex.lab  = 1.9,             # **dimensione etichette assi (xlab/ylab/zlab)**
  cex.axis = 1.4,             # **dimensione valori sui tick**
  cex.main = 1.4              # se usi 'main=' e vuoi più grande
)

## --- LEGEND DISCRETA, attaccata al bordo del grafico ---
par(xpd = NA)  # consenti di disegnare anche sul margine
# Posizione molto vicina al bordo destro/alto del pannello:
lx <- grconvertX(0.75, from = "ndc", to = "user")  # più alto → più vicino al bordo destro
ly <- grconvertY(0.75,  from = "ndc", to = "user")  # regola l'altezza verticale
legend(x = lx, y = ly,
       legend = lvl,
       pch = 15, pt.cex = 2,
       col = unname(cols),
       bty = "n",
       xjust = 1, yjust = 1)

par(op)

ppp_appl_terr_stp <- stpm(data.frame(x = ppp_appl_terr$long, y = ppp_appl_terr$lat, t = ppp_appl_terr$z, m1 = ppp_appl_terr$magn1))
plot(ppp_appl_terr_stp)

ppp_appl_terr_stp$df$m1 <- cut(ppp_appl_terr_stp$df$m1, breaks = c(3, 3.2, 4, 7))
levels(ppp_appl_terr_stp$df$m1) <- c("Mag_low", "Mag_med", "Mag_high")
ppp_appl_terr_stp$df$m1 <- as.factor(ppp_appl_terr_stp$df$m1)

faglie_std <- faglie
faglie_std <- faglie_std[which(faglie_std$x >= min(ppp_appl_terr_stp$df$x) & 
                                 faglie_std$x <= max(ppp_appl_terr_stp$df$x) &
                                 faglie_std$y >= min(ppp_appl_terr_stp$df$y) & 
                                 faglie_std$y <= max(ppp_appl_terr_stp$df$y)),]
faglie_std$x <- (faglie_std$x - min(ppp_appl_terr_stp$df$x))/(max(ppp_appl_terr_stp$df$x) - min(ppp_appl_terr_stp$df$x))
faglie_std$y <- (faglie_std$y - min(ppp_appl_terr_stp$df$y))/(max(ppp_appl_terr_stp$df$y) - min(ppp_appl_terr_stp$df$y))

vul_mod_std <- vul_mod
vul_mod_std <- vul_mod_std[which(vul_mod_std$Longitude >= min(ppp_appl_terr_stp$df$x) & 
                                   vul_mod_std$Longitude <= max(ppp_appl_terr_stp$df$x) &
                                   vul_mod_std$Latitude >= min(ppp_appl_terr_stp$df$y) & 
                                   vul_mod_std$Latitude <= max(ppp_appl_terr_stp$df$y)),]

vul_mod_std$Longitude <- (vul_mod_std$Longitude - min(ppp_appl_terr_stp$df$x))/(max(ppp_appl_terr_stp$df$x) - min(ppp_appl_terr_stp$df$x))
vul_mod_std$Latitude <- (vul_mod_std$Latitude - min(ppp_appl_terr_stp$df$y))/(max(ppp_appl_terr_stp$df$y) - min(ppp_appl_terr_stp$df$y))

ppp_appl_terr_stp$df$x <- (ppp_appl_terr_stp$df$x - min(ppp_appl_terr_stp$df$x))/(max(ppp_appl_terr_stp$df$x) - min(ppp_appl_terr_stp$df$x))
ppp_appl_terr_stp$df$y <- (ppp_appl_terr_stp$df$y - min(ppp_appl_terr_stp$df$y))/(max(ppp_appl_terr_stp$df$y) - min(ppp_appl_terr_stp$df$y))
ppp_appl_terr_stp$df$t <- (ppp_appl_terr_stp$df$t - min(ppp_appl_terr_stp$df$t))/(max(ppp_appl_terr_stp$df$t)-min(ppp_appl_terr_stp$df$t))

terr_model_v1 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = TRUE,
                                     marked = T, spatial.cov = T, mult = 4, ncube = 11, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 3)
terr_model_v11 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = TRUE,
                                      marked = T, spatial.cov = T, mult = 1, ncube = 1, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 3)
terr_model_v2 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = F,
                                     marked = T, spatial.cov = T, mult = 3, ncube = 1, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 2)
terr_model_v3 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = T,
                                     marked = T, spatial.cov = T, mult = 8, ncube = 15, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 2)
terr_model_v4 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = T,
                                     marked = T, spatial.cov = T, mult = 8, ncube = 1, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 2)
terr_model_v5 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = F,
                                     marked = T, spatial.cov = T, mult = 6, ncube = 7, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 2)
terr_model_v6 <- stppm.terremoti.WPI(ppp_appl_terr_stp, formula = ~ x + y + t + s(m1, bs = "re") + Fau, grid = T,
                                     marked = T, spatial.cov = T, mult = 3, ncube = 9, covs = list(Fau = faglie_std), dist.cov = TRUE, seed = 2)


### MOD V1 ###

mod = terr_model_v1
nsim = 200

sim_l1_terr_v1 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = nsim, seed = 1)
sim_l2_terr_v1 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = nsim, seed = 2)
sim_l3_terr_v1 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] - mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = nsim, seed = 3)

sim_terr_ppp_v1 <- setNames(replicate(nsim, data.frame()), paste0("X", 1:nsim))

for(i in 1:nsim) {
  sim_terr_ppp_v1[[i]] <- data.frame(x = c(sim_l1_terr_v1[[i]]$x, sim_l2_terr_v1[[i]]$x, sim_l3_terr_v1[[i]]$x),
                                     y = c(sim_l1_terr_v1[[i]]$y, sim_l2_terr_v1[[i]]$y, sim_l3_terr_v1[[i]]$y),
                                     t = c(sim_l1_terr_v1[[i]]$t, sim_l2_terr_v1[[i]]$t, sim_l3_terr_v1[[i]]$t),
                                     Fau = c(sim_l1_terr_v1[[i]]$Fau, sim_l2_terr_v1[[i]]$Fau, sim_l3_terr_v1[[i]]$Fau),
                                     m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v1[[i]])),
                                                      rep("Mag_med", nrow(sim_l2_terr_v1[[i]])),
                                                      rep("Mag_high", nrow(sim_l3_terr_v1[[i]])))))
}



Khat_data_terr_v1 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v1 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v1$dist), nrow = 100))
for(j in 1:length(sim_terr_ppp_v1)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v1[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v1[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v1$dist)
  Khat_sim_terr_v1[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v1 <- create_curve_set(list(r = Khat_data_terr_v1$dist, obs = Khat_data_terr_v1$Khat, sim_m = t(Khat_sim_terr_v1)))
global_env_terr_v1 <- global_envelope_test(curve_set_terr_v1, alternative = "two.sided")
plot(global_env_terr_v1)

p_v1 <- plot(global_env_terr_v1)
p_v1 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)



### MOD V11 ###

mod = terr_model_v11
nsim = 200

sim_l1_terr_v11 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                       covs = list(Fau = mod$newdata$Fau),
                                       covs_values = list(Fau = faglie_std), nsim = nsim, seed = 1)
sim_l2_terr_v11 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                       covs = list(Fau = mod$newdata$Fau),
                                       covs_values = list(Fau = faglie_std), nsim = nsim, seed = 2)
sim_l3_terr_v11 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                       par = c(mod$IntCoefs[1] - mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                       covs = list(Fau = mod$newdata$Fau),
                                       covs_values = list(Fau = faglie_std), nsim = nsim, seed = 3)

sim_terr_ppp_v11 <- setNames(replicate(nsim, data.frame()), paste0("X", 1:nsim))

for(i in 1:nsim) {
  sim_terr_ppp_v11[[i]] <- data.frame(x = c(sim_l1_terr_v11[[i]]$x, sim_l2_terr_v11[[i]]$x, sim_l3_terr_v11[[i]]$x),
                                      y = c(sim_l1_terr_v11[[i]]$y, sim_l2_terr_v11[[i]]$y, sim_l3_terr_v11[[i]]$y),
                                      t = c(sim_l1_terr_v11[[i]]$t, sim_l2_terr_v11[[i]]$t, sim_l3_terr_v11[[i]]$t),
                                      Fau = c(sim_l1_terr_v11[[i]]$Fau, sim_l2_terr_v11[[i]]$Fau, sim_l3_terr_v11[[i]]$Fau),
                                      m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v11[[i]])),
                                                       rep("Mag_med", nrow(sim_l2_terr_v11[[i]])),
                                                       rep("Mag_high", nrow(sim_l3_terr_v11[[i]])))))
}



Khat_data_terr_v11 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v11 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v11$dist), nrow = 100))
for(j in 1:length(sim_terr_ppp_v11)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v11[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v11[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v11$dist)
  Khat_sim_terr_v11[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v11 <- create_curve_set(list(r = Khat_data_terr_v11$dist, obs = Khat_data_terr_v11$Khat, sim_m = t(Khat_sim_terr_v11)))
global_env_terr_v11 <- global_envelope_test(curve_set_terr_v11, alternative = "two.sided")
plot(global_env_terr_v11)

p_v11 <- plot(global_env_terr_v11)
p_v11 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)




### MOD V2 ###

mod = terr_model_v2

sim_l1_terr_v2 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 1)
sim_l2_terr_v2 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 2)
sim_l3_terr_v2 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 3)

sim_terr_ppp_v2 <- setNames(replicate(200, data.frame()), paste0("X", 1:100))

for(i in 1:200) {
  sim_terr_ppp_v2[[i]] <- data.frame(x = c(sim_l1_terr_v2[[i]]$x, sim_l2_terr_v2[[i]]$x, sim_l3_terr_v2[[i]]$x),
                                     y = c(sim_l1_terr_v2[[i]]$y, sim_l2_terr_v2[[i]]$y, sim_l3_terr_v2[[i]]$y),
                                     t = c(sim_l1_terr_v2[[i]]$t, sim_l2_terr_v2[[i]]$t, sim_l3_terr_v2[[i]]$t),
                                     Fau = c(sim_l1_terr_v2[[i]]$Fau, sim_l2_terr_v2[[i]]$Fau, sim_l3_terr_v2[[i]]$Fau),
                                     m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v2[[i]])),
                                                      rep("Mag_med", nrow(sim_l2_terr_v2[[i]])),
                                                      rep("Mag_high", nrow(sim_l3_terr_v2[[i]])))))
}



Khat_data_terr_v2 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v2 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v2$dist), nrow = 200))
for(j in 1:length(sim_terr_ppp)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v2[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v2[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v2$dist)
  Khat_sim_terr_v2[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v2 <- create_curve_set(list(r = Khat_data_terr_v2$dist, obs = Khat_data_terr_v2$Khat, sim_m = t(Khat_sim_terr_v2)))
global_env_terr_v2 <- global_envelope_test(curve_set_terr_v2, alternative = "two.sided")
plot(global_env_terr_v2)

p_v2 <- plot(global_env_terr_v2)
p_v2 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)


### MOD v3 ###

mod = terr_model_v3

sim_l1_terr_v3 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 1)
sim_l2_terr_v3 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 2)
sim_l3_terr_v3 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 3)

sim_terr_ppp_v3 <- setNames(replicate(200, data.frame()), paste0("X", 1:100))

for(i in 1:200) {
  sim_terr_ppp_v3[[i]] <- data.frame(x = c(sim_l1_terr_v3[[i]]$x, sim_l2_terr_v3[[i]]$x, sim_l3_terr_v3[[i]]$x),
                                     y = c(sim_l1_terr_v3[[i]]$y, sim_l2_terr_v3[[i]]$y, sim_l3_terr_v3[[i]]$y),
                                     t = c(sim_l1_terr_v3[[i]]$t, sim_l2_terr_v3[[i]]$t, sim_l3_terr_v3[[i]]$t),
                                     Fau = c(sim_l1_terr_v3[[i]]$Fau, sim_l2_terr_v3[[i]]$Fau, sim_l3_terr_v3[[i]]$Fau),
                                     m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v3[[i]])),
                                                      rep("Mag_med", nrow(sim_l2_terr_v3[[i]])),
                                                      rep("Mag_high", nrow(sim_l3_terr_v3[[i]])))))
}



Khat_data_terr_v3 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v3 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v3$dist), nrow = 200))
for(j in 1:length(sim_terr_ppp)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v3[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v3[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v3$dist)
  Khat_sim_terr_v3[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v3 <- create_curve_set(list(r = Khat_data_terr_v3$dist, obs = Khat_data_terr_v3$Khat, sim_m = t(Khat_sim_terr_v3)))
global_env_terr_v3 <- global_envelope_test(curve_set_terr_v3, alternative = "two.sided")
plot(global_env_terr_v3)

p_v3 <- plot(global_env_terr_v3)
p_v3 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)



### MOD v4 ###

mod = terr_model_v4

sim_l1_terr_v4 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 1)
sim_l2_terr_v4 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 2)
sim_l3_terr_v4 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 3)

sim_terr_ppp_v4 <- setNames(replicate(200, data.frame()), paste0("X", 1:100))

for(i in 1:200) {
  sim_terr_ppp_v4[[i]] <- data.frame(x = c(sim_l1_terr_v4[[i]]$x, sim_l2_terr_v4[[i]]$x, sim_l3_terr_v4[[i]]$x),
                                     y = c(sim_l1_terr_v4[[i]]$y, sim_l2_terr_v4[[i]]$y, sim_l3_terr_v4[[i]]$y),
                                     t = c(sim_l1_terr_v4[[i]]$t, sim_l2_terr_v4[[i]]$t, sim_l3_terr_v4[[i]]$t),
                                     Fau = c(sim_l1_terr_v4[[i]]$Fau, sim_l2_terr_v4[[i]]$Fau, sim_l3_terr_v4[[i]]$Fau),
                                     m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v4[[i]])),
                                                      rep("Mag_med", nrow(sim_l2_terr_v4[[i]])),
                                                      rep("Mag_high", nrow(sim_l3_terr_v4[[i]])))))
}



Khat_data_terr_v4 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v4 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v4$dist), nrow = 200))
for(j in 1:length(sim_terr_ppp)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v4[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v4[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v4$dist)
  Khat_sim_terr_v4[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v4 <- create_curve_set(list(r = Khat_data_terr_v4$dist, obs = Khat_data_terr_v4$Khat, sim_m = t(Khat_sim_terr_v4)))
global_env_terr_v4 <- global_envelope_test(curve_set_terr_v4, alternative = "two.sided")
plot(global_env_terr_v4)

p_v4 <- plot(global_env_terr_v4)
p_v4 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)



### MOD v5 ###

mod = terr_model_v5

sim_l1_terr_v5 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 1)
sim_l2_terr_v5 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 2)
sim_l3_terr_v5 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 3)

sim_terr_ppp_v5 <- setNames(replicate(200, data.frame()), paste0("X", 1:100))

for(i in 1:200) {
  sim_terr_ppp_v5[[i]] <- data.frame(x = c(sim_l1_terr_v5[[i]]$x, sim_l2_terr_v5[[i]]$x, sim_l3_terr_v5[[i]]$x),
                                     y = c(sim_l1_terr_v5[[i]]$y, sim_l2_terr_v5[[i]]$y, sim_l3_terr_v5[[i]]$y),
                                     t = c(sim_l1_terr_v5[[i]]$t, sim_l2_terr_v5[[i]]$t, sim_l3_terr_v5[[i]]$t),
                                     Fau = c(sim_l1_terr_v5[[i]]$Fau, sim_l2_terr_v5[[i]]$Fau, sim_l3_terr_v5[[i]]$Fau),
                                     m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v5[[i]])),
                                                      rep("Mag_med", nrow(sim_l2_terr_v5[[i]])),
                                                      rep("Mag_high", nrow(sim_l3_terr_v5[[i]])))))
}



Khat_data_terr_v5 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v5 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v5$dist), nrow = 200))
for(j in 1:length(sim_terr_ppp)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v5[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v5[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v5$dist)
  Khat_sim_terr_v5[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v5 <- create_curve_set(list(r = Khat_data_terr_v5$dist, obs = Khat_data_terr_v5$Khat, sim_m = t(Khat_sim_terr_v5)))
global_env_terr_v5 <- global_envelope_test(curve_set_terr_v5, alternative = "two.sided")
plot(global_env_terr_v5)

p_v5 <- plot(global_env_terr_v5)
p_v5 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)



### MOD v6 ###

mod = terr_model_v6

sim_l1_terr_v6 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[6], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 1)
sim_l2_terr_v6 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[7], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 2)
sim_l3_terr_v6 <- rstpp_cov_appl_terr(lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*x + a[3]*y + a[4]*t + a[5]*cov[1])}, 
                                      par = c(mod$IntCoefs[1] + mod$IntCoefs[8], mod$IntCoefs[2], mod$IntCoefs[3], mod$IntCoefs[4], mod$IntCoefs[5]), 
                                      covs = list(Fau = mod$newdata$Fau),
                                      covs_values = list(Fau = faglie_std), nsim = 200, seed = 3)

sim_terr_ppp_v6 <- setNames(replicate(200, data.frame()), paste0("X", 1:100))

for(i in 1:200) {
  sim_terr_ppp_v6[[i]] <- data.frame(x = c(sim_l1_terr_v6[[i]]$x, sim_l2_terr_v6[[i]]$x, sim_l3_terr_v6[[i]]$x),
                                     y = c(sim_l1_terr_v6[[i]]$y, sim_l2_terr_v6[[i]]$y, sim_l3_terr_v6[[i]]$y),
                                     t = c(sim_l1_terr_v6[[i]]$t, sim_l2_terr_v6[[i]]$t, sim_l3_terr_v6[[i]]$t),
                                     Fau = c(sim_l1_terr_v6[[i]]$Fau, sim_l2_terr_v6[[i]]$Fau, sim_l3_terr_v6[[i]]$Fau),
                                     m1 = as.factor(c(rep("Mag_low", nrow(sim_l1_terr_v6[[i]])),
                                                      rep("Mag_med", nrow(sim_l2_terr_v6[[i]])),
                                                      rep("Mag_high", nrow(sim_l3_terr_v6[[i]])))))
}



Khat_data_terr_v6 <- Khat_spatial.3D(stp(ppp_appl_terr_stp$df[,c(1,2,3)]), lambda = mod$l, correction = "translate")

Khat_sim_terr_v6 <- as.data.frame(matrix(NA, ncol = length(Khat_data_terr_v6$dist), nrow = 200))
for(j in 1:length(sim_terr_ppp)){
  pred_int <- exp(predict(mod$mod_global, newdata = sim_terr_ppp_v6[[j]]))
  Khat_sim_terr_j <- Khat_spatial.3D(stp(sim_terr_ppp_v6[[j]][,c(1,2,3)]), lambda = pred_int, correction = "translate", dist = Khat_data_terr_v6$dist)
  Khat_sim_terr_v6[j, ] <- Khat_sim_terr_j$Khat
}

curve_set_terr_v6 <- create_curve_set(list(r = Khat_data_terr_v6$dist, obs = Khat_data_terr_v6$Khat, sim_m = t(Khat_sim_terr_v6)))
global_env_terr_v6 <- global_envelope_test(curve_set_terr_v6, alternative = "two.sided")
plot(global_env_terr_v6)

p_v6 <- plot(global_env_terr_v6)
p_v6 + ggplot2::theme(
  text = ggplot2::element_text(size = 18),
  axis.title = ggplot2::element_text(size = 18),
  axis.text  = ggplot2::element_text(size = 16),
  legend.title = ggplot2::element_text(size = 16),
  legend.text  = ggplot2::element_text(size = 16)
)



### PLOT DEI RISULTATI DEI MODELLI ###

my_cols <- c("#132B43", "#56B1F7", "#8FE388", "#FDE725")  # <--- CAMBIA QUI
pal_fun <- grDevices::colorRampPalette(my_cols)
cols    <- pal_fun(150)

par(mar = c(1, 1, 1, 20))
common_range <- range(c(terr_model_v1$l, terr_model_v11$l, terr_model_v2$l, terr_model_v4$l, terr_model_v5$l, terr_model_v6$l))  # oppure unisci v1 e v2 se vuoi range globale

plot3D::scatter3D(quakes_cat$lon, quakes_cat$lat, quakes_cat$depth,
                  theta = -45, phi = 35,
                  colvar = terr_model_v2$l,
                  col = grDevices::hcl.colors(100, "Inferno", rev = TRUE),
                  clim = common_range,
                  ticktype = "detailed", pch = 20,
                  xlab = "x", ylab = "y", zlab = "z",  cex.lab  = 1.7, cex.axis = 1.4,
                  colkey = list(length = 0.45, width = 0.5, cex.clab = 1.4,  cex.axis = 1.4))

plot3D::scatter3D(quakes_cat$lon, quakes_cat$lat, quakes_cat$depth,
                  theta = -45, phi = 35,
                  colvar = terr_model_v4$l,
                  col = grDevices::hcl.colors(100, "Inferno", rev = TRUE),
                  clim = common_range,
                  ticktype = "detailed", pch = 20,
                  xlab = "x", ylab = "y", zlab = "z",  cex.lab  = 1.7, cex.axis = 1.4,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 1.4,  cex.axis = 1.4))

plot3D::scatter3D(quakes_cat$lon, quakes_cat$lat, quakes_cat$depth,
                  theta = -45, phi = 35,
                  colvar = terr_model_v5$l,
                  col = grDevices::hcl.colors(100, "Inferno", rev = TRUE),
                  clim = common_range,
                  ticktype = "detailed", pch = 20,
                  xlab = "x", ylab = "y", zlab = "z",  cex.lab  = 1.7, cex.axis = 1.4,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 1.4,  cex.axis = 1.4))

plot3D::scatter3D(quakes_cat$lon, quakes_cat$lat, quakes_cat$depth,
                  theta = -45, phi = 35,
                  colvar = terr_model_v6$l,
                  col = grDevices::hcl.colors(100, "Inferno", rev = TRUE),
                  clim = common_range,
                  ticktype = "detailed", pch = 20,
                  xlab = "x", ylab = "y", zlab = "z",  cex.lab  = 1.7, cex.axis = 1.4,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 1.4,  cex.axis = 1.4))

plot3D::scatter3D(quakes_cat$lon, quakes_cat$lat, quakes_cat$depth,
                  theta = -45, phi = 35,
                  colvar = terr_model_v1$l,
                  col = grDevices::hcl.colors(100, "Inferno", rev = TRUE),
                  clim = common_range,
                  ticktype = "detailed", pch = 20,
                  xlab = "x", ylab = "y", zlab = "z",  cex.lab  = 1.7, cex.axis = 1.4,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 1.4,  cex.axis = 1.4))

plot3D::scatter3D(quakes_cat$lon, quakes_cat$lat, quakes_cat$depth,
                  theta = -45, phi = 35,
                  colvar = terr_model_v11$l,
                  col = grDevices::hcl.colors(100, "Inferno", rev = TRUE),
                  clim = common_range,
                  ticktype = "detailed", pch = 20,
                  xlab = "x", ylab = "y", zlab = "z",  cex.lab  = 1.7, cex.axis = 1.4,
                  colkey = list(length = 0.5, width = 0.5, cex.clab = 1.4,  cex.axis = 1.4))


#### FUNCTION FOR THE APPLICATION ####

rstpp_cov_appl_terr <- function(lambda = 500, nsim = 1, seed = 2, verbose = TRUE, 
                                par = NULL, minX = 0, maxX = 1, minY = 0, maxY = 1,
                                minT = 0, maxT = 1, parallel = FALSE, cl = 2, covs, covs_values) 
{
  if (is.numeric(lambda)) {
    par <- log(lambda)
    lambda <- function(x, y, t, a) {
      exp(a[1])
    }
  }
  if (nsim != 1) {
    pp0 <- vector("list", length = nsim)
  }
  set.seed(seed)
  for (i in 1:nsim) {
    if (verbose == T) 
      progressreport(i, nsim)
    
    valX <- ifelse(par[2] > 0, maxX, minX)
    valY <- ifelse(par[3] > 0, maxY, minY)
    valT <- ifelse(par[4] > 0, maxT, minT)
    val_cov1 <- ifelse(par[5] > 0, max(covs[[1]]), min(covs[[1]]))
    lmax <- lambda(valX, valY, valT, cov = c(val_cov1), a = par)
    
    candn <- rpois(1, lmax)
    candx <- runif(candn, min = minX, max = maxX)
    candy <- runif(candn, min = minY, max = maxY)
    candt <- runif(candn, min = minT, max = maxT)
    d <- runif(candn)
    points_sim <- data.frame(x = candx, y = candy, t = candt)
    
    covs_interp <- as.data.frame(matrix(NA, ncol = length(covs), nrow = candn))
    for (j in 1:length(covs)) {
      cov.dist <- compute.dist.cov(points_sim, covs_values[[j]])
      covs_interp[,j] <- cov.dist
    }
    
    lam2 <- numeric(candn)
    for(k in 1:candn) {
      lam2[k] <- as.numeric(lambda(x = candx[k], y = candy[k], t = candt[k], cov = covs_interp[k,], a = par))
    }
    
    keep <- (d < lam2 / lmax)
    n <- sum(keep)
    t1 <- candt[keep]
    t <- t1[order(t1)]
    lon <- candx[keep][order(t1)]
    lat <- candy[keep][order(t1)]
    fau <- covs_interp[keep,1][order(t1)]
    if (nsim != 1) {
      pp0[[i]] <- data.frame(x = lon, y = lat, t = t, Fau = fau)
    }
    else {
      pp0 <- data.frame(x = lon, y = lat, t = t, Fau = fau)
    }
  }
  return(pp0)
}


prova <- rstpp_cov_generalised_v2(
  lambda = function(x, y, t, cov, a) {exp(a[1] + a[2]*cov[1])},
  par = c(comb$mean_param1, comb$mean_param2),
  nsim = 4,
  covs = list(df_covs$cov1, df_covs$cov2, df_covs$cov3))

rstpp_cov_generalised_v2 <- function(lambda = 500, nsim = 1, seed = 2, verbose = TRUE, 
                                     par = NULL, minX = 0, maxX = 1, minY = 0, maxY = 1,
                                     minT = 0, maxT = 1, parallel = FALSE, cl = 2, covs) 
{
  if (is.numeric(lambda)) {
    par <- log(lambda)
    lambda <- function(x, y, t, a) {
      exp(a[1])
    }
  }
  if (nsim != 1) {
    pp0 <- vector("list", length = nsim)
  }
  set.seed(seed)
  for (ri in 1:nsim) {
    if (verbose == T) 
      progressreport(ri, nsim)
    valX <- ifelse(par[2] > 0, maxX, minX)
    valY <- ifelse(par[3] > 0, maxY, minY)
    valT <- ifelse(par[4] > 0, maxT, minT)
    if(length(par) > 4) {val_cov1 <- ifelse(par[5] > 0, max(covs[[1]]$cov1), min(covs[[1]]$cov1))} else {val_cov1 = 1} 
    if(length(par) > 5) {val_cov2 <- ifelse(par[6] > 0, max(covs[[2]]$cov2), min(covs[[2]]$cov2))} else {val_cov2 = 1}
    if(length(par) > 6) {val_cov3 <- ifelse(par[7] > 0, max(covs[[3]]$cov3), min(covs[[3]]$cov3))} else {val_cov3 = 1}
    lmax <- lambda(valX, valY, valT, cov = c(val_cov1, val_cov2, val_cov3), a = par)

    candn <- rpois(1, lmax)
    candx <- runif(candn, min = minX, max = maxX)
    candy <- runif(candn, min = minY, max = maxY)
    candt <- runif(candn, min = minT, max = maxT)
    d <- runif(candn)
    points_sim <- data.frame(x = candx, y = candy, t = candt)
    
    covs_interp <- as.data.frame(matrix(NA, ncol = length(covs), nrow = candn))
    for (rj in 1:length(covs)) {
      covs_data <- covs[[rj]]
      covs_interp[,rj] <- interp3D_prova3(points = points_sim, covs = covs_data, p = 81, d = 3)
    }
    
    lam2 <- numeric(candn)
    for(rk in 1:candn) {
      lam2[rk] <- as.numeric(lambda(x = candx[rk], y = candy[rk], t = candt[rk], cov = covs_interp[rk,], a = par))
    }
    
    keep <- (d < lam2 / lmax)
    n <- sum(keep)
    t1 <- candt[keep]
    t <- t1
    lon <- candx[keep]
    lat <- candy[keep]
    covs_values <- covs_interp[keep,]
    if (nsim != 1) {
      pp0[[ri]] <- data.frame(x = lon, y = lat, t = t, covs_values)
    }
    else {
      pp0 <- data.frame(x = lon, y = lat, t = t, covs_values)
    }
  }
  return(pp0)
}



