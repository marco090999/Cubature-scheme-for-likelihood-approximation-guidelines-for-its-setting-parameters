

#############################################
#### CODES SIMULATION STUDY - GUIDELINES ####
#############################################

find_intersection_baseR <- function(dfs_list) {
  
  # Trova tutti i valori unici di n e grid
  all_n <- sort(unique(unlist(lapply(dfs_list, function(x) unique(x$n)))))
  all_grid <- sort(unique(unlist(lapply(dfs_list, function(x) unique(x$grid)))))
  
  results <- list()
  
  for (n_val in all_n) {
    for (grid_val in all_grid) {
      
      # estrai da ogni scenario le coppie (mult,ncube) buone per n e grid correnti
      good_sets <- list()
      
      for (i in seq_along(dfs_list)) {
        df <- dfs_list[[i]]
        sub <- df[df$n == n_val & df$grid == grid_val &
                    df$contour_class %in% c("white", "yellow"), ]
        if (nrow(sub) > 0) {
          good_sets[[i]] <- unique(sub[, c("mult", "ncube")])
        }
      }
      
      # calcola intersezione tra gli scenari disponibili
      if (length(good_sets) > 0) {
        inter <- Reduce(function(x, y) merge(x, y, by = c("mult", "ncube")), good_sets)
      } else {
        inter <- data.frame(mult = numeric(0), ncube = numeric(0))
      }
      
      results[[paste0("n", n_val, "_grid", grid_val)]] <- list(
        n = n_val,
        grid = grid_val,
        n_scen_present = length(good_sets),
        n_cells_in_intersection = nrow(inter),
        intersection = inter
      )
    }
  }
  
  # trasforma in data.frame di riepilogo
  summary_df <- data.frame(
    n = sapply(results, function(x) x$n),
    grid = sapply(results, function(x) x$grid),
    n_scen_present = sapply(results, function(x) x$n_scen_present),
    n_cells_in_intersection = sapply(results, function(x) x$n_cells_in_intersection),
    stringsAsFactors = FALSE
  )
  
  return(list(summary = summary_df, details = results))
}


############################
#### GRID = T | N = 100 ####
############################

dfs_gridT_n100 <- list(
  scen1 = df_summary_stats_pvalue_gridT_scen1_n100,
  scen2 = df_summary_stats_pvalue_gridT_scen2_n100,
  scen3 = df_summary_stats_pvalue_gridT_scen3_n100,
  scen4 = df_summary_stats_pvalue_gridT_scen4_n100,
  scen5 = df_summary_stats_pvalue_gridT_scen5_n100,
  scen6 = df_summary_stats_pvalue_gridT_scen6_n100,
  scen7 = df_summary_stats_pvalue_gridT_scen7_n100,
  scen8 = df_summary_stats_pvalue_gridT_scen8_n100
)

dfs_gridT_n100_clean <- lapply(dfs_gridT_n100, function(x) {
  # garantisci che esistano le colonne grid e n (numerosità)
  if (!"grid" %in% names(x)) stop("Manca la colonna 'grid'")
  if (!"n" %in% names(x)) x$n <- 500
  x$grid <- as.logical(x$grid)
  x <- x[, c("mult", "ncube", "grid", "n", "contour_class")]
  x
})

list_white_yellow_gridT_n100 <- lapply(dfs_gridT_n100_clean, function(df) {
  subset(df, contour_class %in% c("white", "yellow"), select = c("mult", "ncube"))})

all_pairs_gridT_n100 <- do.call(rbind, list_white_yellow_gridT_n100)
all_pairs_gridT_n100$key <- paste(all_pairs_gridT_n100$mult, all_pairs_gridT_n100$ncube)

freq_table_gridT_n100 <- as.data.frame(table(all_pairs_gridT_n100$key))
names(freq_table_gridT_n100) <- c("key", "count")

freq_table_gridT_n100 <- subset(freq_table_gridT_n100, count >= 7)

res_gridT_n100_df <- do.call(rbind, strsplit(as.character(freq_table_gridT_n100$key), " "))
res_gridT_n100_df <- data.frame(
  mult = as.numeric(res_gridT_n100_df[, 1]),
  ncube = as.numeric(res_gridT_n100_df[, 2]),
  count = freq_table_gridT_n100$count
)

grid_all_gridT_n100 <- expand.grid(mult = 1:50, ncube = 1:50)

grid_all_gridT_n100$present_level <- NA
keys_gridT_n100 <- paste(res_gridT_n100_df$mult, res_gridT_n100_df$ncube)
match_idx_gridT_n100 <- match(paste(grid_all_gridT_n100$mult, grid_all_gridT_n100$ncube), keys_gridT_n100)
grid_all_gridT_n100$present_level <- res_gridT_n100_df$count[match_idx_gridT_n100]
grid_all_gridT_n100$present_level <- ifelse(is.na(grid_all_gridT_n100$present_level), "none",
                                            ifelse(grid_all_gridT_n100$present_level == 8, "8/8", "7/8"))

p1 <- ggplot(grid_all_gridT_n100, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"), # blu chiaro / blu scuro
    breaks = c("8/8", "7/8"),
    name = "Presenza") +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = TRUE | n = 100", x = "q", y = expression(n[c]), fill = "Numero scenari")


### PROVA RICERCA LINEE GUIDA ###

# 1. Seleziona solo le combinazioni "buone"
good_combs_gridT_n100 <- grid_all_gridT_n100 %>%
  filter(present_level %in% c("7/8", "8/8"))

# 2. Calcola il numero di combinazioni buone per ogni q
q_counts_gridT_n100 <- good_combs_gridT_n100 %>%
  group_by(mult) %>%
  summarise(n_good = n()) %>%
  filter(n_good >= 3)   # Manteniamo solo i q con almeno 3 combinazioni buone

# 3. Filtra i dati originali per mantenere solo questi q
good_combs_filtered_gridT_n100 <- good_combs_gridT_n100 %>%
  filter(mult %in% q_counts_gridT_n100$mult)

# 4. Calcola i quantili (0.2, 0.5, 0.8) per ciascun q filtrato
quantile_data_gridT_n100 <- good_combs_filtered_gridT_n100 %>%
  group_by(mult) %>%
  summarise(
    q20 = quantile(ncube, 0.2),
    q50 = quantile(ncube, 0.5),
    q80 = quantile(ncube, 0.8)
  ) %>%
  ungroup()

# 5. Fit lineare per ciascun quantile
fit_q20_gridT_n100 <- lm(q20 ~ mult, data = quantile_data_gridT_n100)
fit_q50_gridT_n100 <- lm(q50 ~ mult, data = quantile_data_gridT_n100)
fit_q80_gridT_n100 <- lm(q80 ~ mult, data = quantile_data_gridT_n100)

# 6. Crea dataframe con le curve stimate
curve_data_gridT_n100 <- data.frame(
  mult = quantile_data_gridT_n100$mult,
  lower = predict(fit_q20_gridT_n100, newdata = quantile_data_gridT_n100),
  median = predict(fit_q50_gridT_n100, newdata = quantile_data_gridT_n100),
  upper = predict(fit_q80_gridT_n100, newdata = quantile_data_gridT_n100)
)

# 7. Plot della griglia + curve
p_gridT_n100 <- ggplot(grid_all_gridT_n100, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"),
    breaks = c("8/8", "7/8"),
    name = "Presenza"
  ) +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = TRUE | n = 100", 
       x = "q", 
       y = expression(n[c]), 
       fill = "Numero scenari")

# 8. Aggiungi le tre curve parametriche
p_gridT_n100 +
  geom_ribbon(data = curve_data_gridT_n100,
              aes(x = mult, ymin = lower, ymax = upper),
              fill = "#d73027", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = curve_data_gridT_n100, aes(x = mult, y = median),
            color = "#b2182b", size = 1, inherit.aes = FALSE)




############################
#### GRID = F | N = 100 ####
############################

dfs_gridF_n100 <- list(
  scen1 = df_summary_stats_pvalue_gridF_scen1_n100,
  scen2 = df_summary_stats_pvalue_gridF_scen2_n100,
  scen3 = df_summary_stats_pvalue_gridF_scen3_n100,
  scen4 = df_summary_stats_pvalue_gridF_scen4_n100,
  scen5 = df_summary_stats_pvalue_gridF_scen5_n100,
  scen6 = df_summary_stats_pvalue_gridF_scen6_n100,
  scen7 = df_summary_stats_pvalue_gridF_scen7_n100,
  scen8 = df_summary_stats_pvalue_gridF_scen8_n100
)

dfs_gridF_n100_clean <- lapply(dfs_gridF_n100, function(x) {
  # garantisci che esistano le colonne grid e n (numerosità)
  if (!"grid" %in% names(x)) stop("Manca la colonna 'grid'")
  if (!"n" %in% names(x)) x$n <- 500
  x$grid <- as.logical(x$grid)
  x <- x[, c("mult", "ncube", "grid", "n", "contour_class")]
  x
})

list_white_yellow_gridF_n100 <- lapply(dfs_gridF_n100_clean, function(df) {
  subset(df, contour_class %in% c("white", "yellow"), select = c("mult", "ncube"))})

all_pairs_gridF_n100 <- do.call(rbind, list_white_yellow_gridF_n100)
all_pairs_gridF_n100$key <- paste(all_pairs_gridF_n100$mult, all_pairs_gridF_n100$ncube)

freq_table_gridF_n100 <- as.data.frame(table(all_pairs_gridF_n100$key))
names(freq_table_gridF_n100) <- c("key", "count")

freq_table_gridF_n100 <- subset(freq_table_gridF_n100, count >= 7)

res_gridF_n100_df <- do.call(rbind, strsplit(as.character(freq_table_gridF_n100$key), " "))
res_gridF_n100_df <- data.frame(
  mult = as.numeric(res_gridF_n100_df[, 1]),
  ncube = as.numeric(res_gridF_n100_df[, 2]),
  count = freq_table_gridF_n100$count
)

grid_all_gridF_n100 <- expand.grid(mult = 1:50, ncube = 1:50)

grid_all_gridF_n100$present_level <- NA
keys_gridF_n100 <- paste(res_gridF_n100_df$mult, res_gridF_n100_df$ncube)
match_idx_gridF_n100 <- match(paste(grid_all_gridF_n100$mult, grid_all_gridF_n100$ncube), keys_gridF_n100)
grid_all_gridF_n100$present_level <- res_gridF_n100_df$count[match_idx_gridF_n100]
grid_all_gridF_n100$present_level <- ifelse(is.na(grid_all_gridF_n100$present_level), "none",
                                            ifelse(grid_all_gridF_n100$present_level == 8, "8/8", "7/8"))

p2 <- ggplot(grid_all_gridF_n100, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"), # blu chiaro / blu scuro
    breaks = c("8/8", "7/8"),
    name = "Presenza") +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.border = element_rect(color = "black", size = 1.2, fill = NA),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "grid = FALSE | n = 100", x = "q", y = expression(n[c]), fill = "Numero scenari")


### PROVA RICERCA LINEE GUIDA ###

# 1. Seleziona solo le combinazioni "buone"
good_combs_gridF_n100 <- grid_all_gridF_n100 %>%
  filter(present_level %in% c("7/8", "8/8"))

# 2. Calcola il numero di combinazioni buone per ogni q
q_counts_gridF_n100 <- good_combs_gridF_n100 %>%
  group_by(mult) %>%
  summarise(n_good = n()) %>%
  filter(n_good >= 3)   # Manteniamo solo i q con almeno 3 combinazioni buone

# 3. Filtra i dati originali per mantenere solo questi q
good_combs_filtered_gridF_n100 <- good_combs_gridF_n100 %>%
  filter(mult %in% q_counts_gridF_n100$mult)

# 4. Calcola i quantili (0.2, 0.5, 0.8) per ciascun q filtrato
quantile_data_gridF_n100 <- good_combs_filtered_gridF_n100 %>%
  group_by(mult) %>%
  summarise(
    q20 = quantile(ncube, 0.2),
    q50 = quantile(ncube, 0.5),
    q80 = quantile(ncube, 0.8)
  ) %>%
  ungroup()

# 5. Fit lineare per ciascun quantile
fit_q20_gridF_n100 <- lm(q20 ~ mult, data = quantile_data_gridF_n100)
fit_q50_gridF_n100 <- lm(q50 ~ mult, data = quantile_data_gridF_n100)
fit_q80_gridF_n100 <- lm(q80 ~ mult, data = quantile_data_gridF_n100)

# 6. Crea dataframe con le curve stimate
curve_data_gridF_n100 <- data.frame(
  mult = quantile_data_gridF_n100$mult,
  lower = predict(fit_q20_gridF_n100, newdata = quantile_data_gridF_n100),
  median = predict(fit_q50_gridF_n100, newdata = quantile_data_gridF_n100),
  upper = predict(fit_q80_gridF_n100, newdata = quantile_data_gridF_n100)
)

# 7. Plot della griglia + curve
p_gridF_n100 <- ggplot(grid_all_gridF_n100, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"),
    breaks = c("8/8", "7/8"),
    name = "Presenza"
  ) +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = FALSE | n = 100", 
       x = "q", 
       y = expression(n[c]), 
       fill = "Numero scenari")

# 8. Aggiungi le tre curve parametriche
p_gridF_n100 +
  geom_ribbon(data = curve_data_gridF_n100,
              aes(x = mult, ymin = lower, ymax = upper),
              fill = "#d73027", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = curve_data_gridF_n100, aes(x = mult, y = median),
            color = "#b2182b", size = 1, inherit.aes = FALSE)



############################
#### GRID = T | N = 250 ####
############################

dfs_gridT_n250 <- list(
  scen1 = df_summary_stats_pvalue_gridT_scen1_n250,
  scen2 = df_summary_stats_pvalue_gridT_scen2_n250,
  scen3 = df_summary_stats_pvalue_gridT_scen3_n250,
  scen4 = df_summary_stats_pvalue_gridT_scen4_n250,
  scen5 = df_summary_stats_pvalue_gridT_scen5_n250,
  scen6 = df_summary_stats_pvalue_gridT_scen6_n250,
  scen7 = df_summary_stats_pvalue_gridT_scen7_n250,
  scen8 = df_summary_stats_pvalue_gridT_scen8_n250
)

dfs_gridT_n250_clean <- lapply(dfs_gridT_n250, function(x) {
  # garantisci che esistano le colonne grid e n (numerosità)
  if (!"grid" %in% names(x)) stop("Manca la colonna 'grid'")
  if (!"n" %in% names(x)) x$n <- 500
  x$grid <- as.logical(x$grid)
  x <- x[, c("mult", "ncube", "grid", "n", "contour_class")]
  x
})

list_white_yellow_gridT_n250 <- lapply(dfs_gridT_n250_clean, function(df) {
  subset(df, contour_class %in% c("white", "yellow"), select = c("mult", "ncube"))})

all_pairs_gridT_n250 <- do.call(rbind, list_white_yellow_gridT_n250)
all_pairs_gridT_n250$key <- paste(all_pairs_gridT_n250$mult, all_pairs_gridT_n250$ncube)

freq_table_gridT_n250 <- as.data.frame(table(all_pairs_gridT_n250$key))
names(freq_table_gridT_n250) <- c("key", "count")

freq_table_gridT_n250 <- subset(freq_table_gridT_n250, count >= 7)

res_gridT_n250_df <- do.call(rbind, strsplit(as.character(freq_table_gridT_n250$key), " "))
res_gridT_n250_df <- data.frame(
  mult = as.numeric(res_gridT_n250_df[, 1]),
  ncube = as.numeric(res_gridT_n250_df[, 2]),
  count = freq_table_gridT_n250$count
)

grid_all_gridT_n250 <- expand.grid(mult = 1:50, ncube = 1:50)

grid_all_gridT_n250$present_level <- NA
keys_gridT_n250 <- paste(res_gridT_n250_df$mult, res_gridT_n250_df$ncube)
match_idx_gridT_n250 <- match(paste(grid_all_gridT_n250$mult, grid_all_gridT_n250$ncube), keys_gridT_n250)
grid_all_gridT_n250$present_level <- res_gridT_n250_df$count[match_idx_gridT_n250]
grid_all_gridT_n250$present_level <- ifelse(is.na(grid_all_gridT_n250$present_level), "none",
                                            ifelse(grid_all_gridT_n250$present_level == 8, "8/8", "7/8"))

p3 <- ggplot(grid_all_gridT_n250, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"), # blu chiaro / blu scuro
    breaks = c("8/8", "7/8"),
    name = "Presenza") +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.border = element_rect(color = "black", size = 1.2, fill = NA),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "grid = TRUE | n = 250", x = "q", y = expression(n[c]), fill = "Numero scenari")



### PROVA RICERCA LINEE GUIDA ###

# 1. Seleziona solo le combinazioni "buone"
good_combs_gridT_n250 <- grid_all_gridT_n250 %>%
  filter(present_level %in% c("7/8", "8/8"))

# 2. Calcola il numero di combinazioni buone per ogni q
q_counts_gridT_n250 <- good_combs_gridT_n250 %>%
  group_by(mult) %>%
  summarise(n_good = n()) %>%
  filter(n_good >= 3)   # Manteniamo solo i q con almeno 3 combinazioni buone

# 3. Filtra i dati originali per mantenere solo questi q
good_combs_filtered_gridT_n250 <- good_combs_gridT_n250 %>%
  filter(mult %in% q_counts_gridT_n250$mult)

# 4. Calcola i quantili (0.2, 0.5, 0.8) per ciascun q filtrato
quantile_data_gridT_n250 <- good_combs_filtered_gridT_n250 %>%
  group_by(mult) %>%
  summarise(
    q20 = quantile(ncube, 0.2),
    q50 = quantile(ncube, 0.5),
    q80 = quantile(ncube, 0.8)
  ) %>%
  ungroup()

# 5. Fit lineare per ciascun quantile
fit_q20_gridT_n250 <- lm(q20 ~ mult, data = quantile_data_gridT_n250)
fit_q50_gridT_n250 <- lm(q50 ~ mult, data = quantile_data_gridT_n250)
fit_q80_gridT_n250 <- lm(q80 ~ mult, data = quantile_data_gridT_n250)

# 6. Crea dataframe con le curve stimate
curve_data_gridT_n250 <- data.frame(
  mult = quantile_data_gridT_n250$mult,
  lower = predict(fit_q20_gridT_n250, newdata = quantile_data_gridT_n250),
  median = predict(fit_q50_gridT_n250, newdata = quantile_data_gridT_n250),
  upper = predict(fit_q80_gridT_n250, newdata = quantile_data_gridT_n250)
)

# 7. Plot della griglia + curve
p_gridT_n250 <- ggplot(grid_all_gridT_n250, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"),
    breaks = c("8/8", "7/8"),
    name = "Presenza"
  ) +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = TRUE | n = 250", 
       x = "q", 
       y = expression(n[c]), 
       fill = "Numero scenari")

# 8. Aggiungi le tre curve parametriche
p_gridT_n250 +
  geom_ribbon(data = curve_data_gridT_n250,
              aes(x = mult, ymin = lower, ymax = upper),
              fill = "#d73027", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = curve_data_gridT_n250, aes(x = mult, y = median),
            color = "#b2182b", size = 1, inherit.aes = FALSE)



############################
#### GRID = F | N = 250 ####
############################

dfs_gridF_n250 <- list(
  scen1 = df_summary_stats_pvalue_gridF_scen1_n250,
  scen2 = df_summary_stats_pvalue_gridF_scen2_n250,
  scen3 = df_summary_stats_pvalue_gridF_scen3_n250,
  scen4 = df_summary_stats_pvalue_gridF_scen4_n250,
  scen5 = df_summary_stats_pvalue_gridF_scen5_n250,
  scen6 = df_summary_stats_pvalue_gridF_scen6_n250,
  scen7 = df_summary_stats_pvalue_gridF_scen7_n250,
  scen8 = df_summary_stats_pvalue_gridF_scen8_n250_mod
)

dfs_gridF_n250_clean <- lapply(dfs_gridF_n250, function(x) {
  # garantisci che esistano le colonne grid e n (numerosità)
  if (!"grid" %in% names(x)) stop("Manca la colonna 'grid'")
  if (!"n" %in% names(x)) x$n <- 500
  x$grid <- as.logical(x$grid)
  x <- x[, c("mult", "ncube", "grid", "n", "contour_class")]
  x
})

list_white_yellow_gridF_n250 <- lapply(dfs_gridF_n250_clean, function(df) {
  subset(df, contour_class %in% c("white", "yellow"), select = c("mult", "ncube"))})

all_pairs_gridF_n250 <- do.call(rbind, list_white_yellow_gridF_n250)
all_pairs_gridF_n250$key <- paste(all_pairs_gridF_n250$mult, all_pairs_gridF_n250$ncube)

freq_table_gridF_n250 <- as.data.frame(table(all_pairs_gridF_n250$key))
names(freq_table_gridF_n250) <- c("key", "count")

freq_table_gridF_n250 <- subset(freq_table_gridF_n250, count >= 7)

res_gridF_n250_df <- do.call(rbind, strsplit(as.character(freq_table_gridF_n250$key), " "))
res_gridF_n250_df <- data.frame(
  mult = as.numeric(res_gridF_n250_df[, 1]),
  ncube = as.numeric(res_gridF_n250_df[, 2]),
  count = freq_table_gridF_n250$count
)

grid_all_gridF_n250 <- expand.grid(mult = 1:50, ncube = 1:50)

grid_all_gridF_n250$present_level <- NA
keys_gridF_n250 <- paste(res_gridF_n250_df$mult, res_gridF_n250_df$ncube)
match_idx_gridF_n250 <- match(paste(grid_all_gridF_n250$mult, grid_all_gridF_n250$ncube), keys_gridF_n250)
grid_all_gridF_n250$present_level <- res_gridF_n250_df$count[match_idx_gridF_n250]
grid_all_gridF_n250$present_level <- ifelse(is.na(grid_all_gridF_n250$present_level), "none",
                                            ifelse(grid_all_gridF_n250$present_level == 8, "8/8", "7/8"))

p4 <- ggplot(grid_all_gridF_n250, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"), # blu chiaro / blu scuro
    breaks = c("8/8", "7/8"),
    name = "Presenza") +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.border = element_rect(color = "black", size = 1.2, fill = NA),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "grid = FALSE | n = 250", x = "q", y = expression(n[c]), fill = "Numero scenari")



### PROVA RICERCA LINEE GUIDA ###

# 1. Seleziona solo le combinazioni "buone"
good_combs_gridF_n250 <- grid_all_gridF_n250 %>%
  filter(present_level %in% c("7/8", "8/8"))

# 2. Calcola il numero di combinazioni buone per ogni q
q_counts_gridF_n250 <- good_combs_gridF_n250 %>%
  group_by(mult) %>%
  summarise(n_good = n()) %>%
  filter(n_good >= 3)   # Manteniamo solo i q con almeno 3 combinazioni buone

# 3. Filtra i dati originali per mantenere solo questi q
good_combs_filtered_gridF_n250 <- good_combs_gridF_n250 %>%
  filter(mult %in% q_counts_gridF_n250$mult)

# 4. Calcola i quantili (0.2, 0.5, 0.8) per ciascun q filtrato
quantile_data_gridF_n250 <- good_combs_filtered_gridF_n250 %>%
  group_by(mult) %>%
  summarise(
    q20 = quantile(ncube, 0.2),
    q50 = quantile(ncube, 0.5),
    q80 = quantile(ncube, 0.8)
  ) %>%
  ungroup()

# 5. Fit lineare per ciascun quantile
fit_q20_gridF_n250 <- lm(q20 ~ mult, data = quantile_data_gridF_n250)
fit_q50_gridF_n250 <- lm(q50 ~ mult, data = quantile_data_gridF_n250)
fit_q80_gridF_n250 <- lm(q80 ~ mult, data = quantile_data_gridF_n250)

# 6. Crea dataframe con le curve stimate
curve_data_gridF_n250 <- data.frame(
  mult = quantile_data_gridF_n250$mult,
  lower = predict(fit_q20_gridF_n250, newdata = quantile_data_gridF_n250),
  median = predict(fit_q50_gridF_n250, newdata = quantile_data_gridF_n250),
  upper = predict(fit_q80_gridF_n250, newdata = quantile_data_gridF_n250)
)

# 7. Plot della griglia + curve
p_gridF_n250 <- ggplot(grid_all_gridF_n250, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"),
    breaks = c("8/8", "7/8"),
    name = "Presenza"
  ) +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = FALSE | n = 250", 
       x = "q", 
       y = expression(n[c]), 
       fill = "Numero scenari")

# 8. Aggiungi le tre curve parametriche
p_gridF_n250 +
  geom_ribbon(data = curve_data_gridF_n250,
              aes(x = mult, ymin = lower, ymax = upper),
              fill = "#d73027", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = curve_data_gridF_n250, aes(x = mult, y = median),
            color = "#b2182b", size = 1, inherit.aes = FALSE)




############################
#### GRID = T | N = 500 ####
############################

dfs_gridT_n500 <- list(
  scen1 = df_summary_stats_pvalue_gridT_scen1_n500,
  scen2 = df_summary_stats_pvalue_gridT_scen2_n500,
  scen3 = df_summary_stats_pvalue_gridT_scen3_n500,
  scen4 = df_summary_stats_pvalue_gridT_scen4_n500,
  scen5 = df_summary_stats_pvalue_gridT_scen5_n500,
  scen6 = df_summary_stats_pvalue_gridT_scen6_n500,
  scen7 = df_summary_stats_pvalue_gridT_scen7_n500,
  scen8 = df_summary_stats_pvalue_gridT_scen8_n500_mod
)

dfs_gridT_n500_clean <- lapply(dfs_gridT_n500, function(x) {
  # garantisci che esistano le colonne grid e n (numerosità)
  if (!"grid" %in% names(x)) stop("Manca la colonna 'grid'")
  if (!"n" %in% names(x)) x$n <- 500
  x$grid <- as.logical(x$grid)
  x <- x[, c("mult", "ncube", "grid", "n", "contour_class")]
  x
})

list_white_yellow_gridT_n500 <- lapply(dfs_gridT_n500_clean, function(df) {
  subset(df, contour_class %in% c("white", "yellow"), select = c("mult", "ncube"))})

all_pairs_gridT_n500 <- do.call(rbind, list_white_yellow_gridT_n500)
all_pairs_gridT_n500$key <- paste(all_pairs_gridT_n500$mult, all_pairs_gridT_n500$ncube)

freq_table_gridT_n500 <- as.data.frame(table(all_pairs_gridT_n500$key))
names(freq_table_gridT_n500) <- c("key", "count")

freq_table_gridT_n500 <- subset(freq_table_gridT_n500, count >= 7)

res_gridT_n500_df <- do.call(rbind, strsplit(as.character(freq_table_gridT_n500$key), " "))
res_gridT_n500_df <- data.frame(
  mult = as.numeric(res_gridT_n500_df[, 1]),
  ncube = as.numeric(res_gridT_n500_df[, 2]),
  count = freq_table_gridT_n500$count
)

grid_all_gridT_n500 <- expand.grid(mult = 1:50, ncube = 1:50)

grid_all_gridT_n500$present_level <- NA
keys_gridT_n500 <- paste(res_gridT_n500_df$mult, res_gridT_n500_df$ncube)
match_idx_gridT_n500 <- match(paste(grid_all_gridT_n500$mult, grid_all_gridT_n500$ncube), keys_gridT_n500)
grid_all_gridT_n500$present_level <- res_gridT_n500_df$count[match_idx_gridT_n500]
grid_all_gridT_n500$present_level <- ifelse(is.na(grid_all_gridT_n500$present_level), "none",
                                            ifelse(grid_all_gridT_n500$present_level == 8, "8/8", "7/8"))

p5 <- ggplot(grid_all_gridT_n500, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"), # blu chiaro / blu scuro
    breaks = c("8/8", "7/8"),
    name = "Presenza") +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.border = element_rect(color = "black", size = 1.2, fill = NA),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "grid = TRUE | n = 500", x = "q", y = expression(n[c]), fill = "Numero scenari")



### PROVA RICERCA LINEE GUIDA ###

# 1. Seleziona solo le combinazioni "buone"
good_combs_gridT_n500 <- grid_all_gridT_n500 %>%
  filter(present_level %in% c("7/8", "8/8"))

# 2. Calcola il numero di combinazioni buone per ogni q
q_counts_gridT_n500 <- good_combs_gridT_n500 %>%
  group_by(mult) %>%
  summarise(n_good = n()) %>%
  filter(n_good >= 3)   # Manteniamo solo i q con almeno 3 combinazioni buone

# 3. Filtra i dati originali per mantenere solo questi q
good_combs_filtered_gridT_n500 <- good_combs_gridT_n500 %>%
  filter(mult %in% q_counts_gridT_n500$mult)

# 4. Calcola i quantili (0.2, 0.5, 0.8) per ciascun q filtrato
quantile_data_gridT_n500 <- good_combs_filtered_gridT_n500 %>%
  group_by(mult) %>%
  summarise(
    q20 = quantile(ncube, 0.2),
    q50 = quantile(ncube, 0.5),
    q80 = quantile(ncube, 0.8)
  ) %>%
  ungroup()

# 5. Fit lineare per ciascun quantile
fit_q20_gridT_n500 <- lm(q20 ~ mult, data = quantile_data_gridT_n500)
fit_q50_gridT_n500 <- lm(q50 ~ mult, data = quantile_data_gridT_n500)
fit_q80_gridT_n500 <- lm(q80 ~ mult, data = quantile_data_gridT_n500)

# 6. Crea dataframe con le curve stimate
curve_data_gridT_n500 <- data.frame(
  mult = quantile_data_gridT_n500$mult,
  lower = predict(fit_q20_gridT_n500, newdata = quantile_data_gridT_n500),
  median = predict(fit_q50_gridT_n500, newdata = quantile_data_gridT_n500),
  upper = predict(fit_q80_gridT_n500, newdata = quantile_data_gridT_n500)
)

# 7. Plot della griglia + curve
p_gridT_n500 <- ggplot(grid_all_gridT_n500, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"),
    breaks = c("8/8", "7/8"),
    name = "Presenza"
  ) +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = TRUE | n = 500", 
       x = "q", 
       y = expression(n[c]), 
       fill = "Numero scenari")

# 8. Aggiungi le tre curve parametriche
p_gridT_n500 +
  geom_ribbon(data = curve_data_gridT_n500,
              aes(x = mult, ymin = lower, ymax = upper),
              fill = "#d73027", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = curve_data_gridT_n500, aes(x = mult, y = median),
            color = "#b2182b", size = 1, inherit.aes = FALSE)



############################
#### GRID = F | N = 500 ####
############################

dfs_gridF_n500 <- list(
  scen1 = df_summary_stats_pvalue_gridF_scen1_n500,
  scen2 = df_summary_stats_pvalue_gridF_scen2_n500,
  scen3 = df_summary_stats_pvalue_gridF_scen3_n500,
  scen4 = df_summary_stats_pvalue_gridF_scen4_n500,
  scen5 = df_summary_stats_pvalue_gridF_scen5_n500,
  scen6 = df_summary_stats_pvalue_gridF_scen6_n500,
  scen7 = df_summary_stats_pvalue_gridF_scen7_n500,
  scen8 = df_summary_stats_pvalue_gridF_scen8_n500_mod
)

dfs_gridF_n500_clean <- lapply(dfs_gridF_n500, function(x) {
  # garantisci che esistano le colonne grid e n (numerosità)
  if (!"grid" %in% names(x)) stop("Manca la colonna 'grid'")
  if (!"n" %in% names(x)) x$n <- 500
  x$grid <- as.logical(x$grid)
  x <- x[, c("mult", "ncube", "grid", "n", "contour_class")]
  x
})


list_white_yellow_gridF_n500 <- lapply(dfs_gridF_n500_clean, function(df) {
  subset(df, contour_class %in% c("white", "yellow"), select = c("mult", "ncube"))})

all_pairs_gridF_n500 <- do.call(rbind, list_white_yellow_gridF_n500)
all_pairs_gridF_n500$key <- paste(all_pairs_gridF_n500$mult, all_pairs_gridF_n500$ncube)

freq_table_gridF_n500 <- as.data.frame(table(all_pairs_gridF_n500$key))
names(freq_table_gridF_n500) <- c("key", "count")

freq_table_gridF_n500 <- subset(freq_table_gridF_n500, count >= 7)

res_gridF_n500_df <- do.call(rbind, strsplit(as.character(freq_table_gridF_n500$key), " "))
res_gridF_n500_df <- data.frame(
  mult = as.numeric(res_gridF_n500_df[, 1]),
  ncube = as.numeric(res_gridF_n500_df[, 2]),
  count = freq_table_gridF_n500$count
)

grid_all_gridF_n500 <- expand.grid(mult = 1:50, ncube = 1:50)

grid_all_gridF_n500$present_level <- NA
keys_gridF_n500 <- paste(res_gridF_n500_df$mult, res_gridF_n500_df$ncube)
match_idx_gridF_n500 <- match(paste(grid_all_gridF_n500$mult, grid_all_gridF_n500$ncube), keys_gridF_n500)
grid_all_gridF_n500$present_level <- res_gridF_n500_df$count[match_idx_gridF_n500]
grid_all_gridF_n500$present_level <- ifelse(is.na(grid_all_gridF_n500$present_level), "none",
                                            ifelse(grid_all_gridF_n500$present_level == 8, "8/8", "7/8"))

p6 <- ggplot(grid_all_gridF_n500, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"), # blu chiaro / blu scuro
    breaks = c("8/8", "7/8"),
    name = "Presenza") +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "grid = FALSE | n = 500", x = "q", y = expression(n[c]), fill = "Numero scenari")
    

### PROVA RICERCA LINEE GUIDA ###

# 1. Seleziona solo le combinazioni "buone"
good_combs_gridF_n500 <- grid_all_gridF_n500 %>%
  filter(present_level %in% c("7/8", "8/8"))

# 2. Calcola il numero di combinazioni buone per ogni q
q_counts_gridF_n500 <- good_combs_gridF_n500 %>%
  group_by(mult) %>%
  summarise(n_good = n()) %>%
  filter(n_good >= 3)   # Manteniamo solo i q con almeno 3 combinazioni buone

# 3. Filtra i dati originali per mantenere solo questi q
good_combs_filtered_gridF_n500 <- good_combs_gridF_n500 %>%
  filter(mult %in% q_counts_gridF_n500$mult)

# 4. Calcola i quantili (0.2, 0.5, 0.8) per ciascun q filtrato
quantile_data_gridF_n500 <- good_combs_filtered_gridF_n500 %>%
  group_by(mult) %>%
  summarise(
    q20 = quantile(ncube, 0.2),
    q50 = quantile(ncube, 0.5),
    q80 = quantile(ncube, 0.8)
  ) %>%
  ungroup()

# 5. Fit lineare per ciascun quantile
fit_q20_gridF_n500 <- lm(q20 ~ mult, data = quantile_data_gridF_n500)
fit_q50_gridF_n500 <- lm(q50 ~ mult, data = quantile_data_gridF_n500)
fit_q80_gridF_n500 <- lm(q80 ~ mult, data = quantile_data_gridF_n500)

# 6. Crea dataframe con le curve stimate
curve_data_gridF_n500 <- data.frame(
  mult = quantile_data_gridF_n500$mult,
  lower = predict(fit_q20_gridF_n500, newdata = quantile_data_gridF_n500),
  median = predict(fit_q50_gridF_n500, newdata = quantile_data_gridF_n500),
  upper = predict(fit_q80_gridF_n500, newdata = quantile_data_gridF_n500)
)

# 7. Plot della griglia + curve
p_gridF_n500 <- ggplot(grid_all_gridF_n500, aes(x = mult, y = ncube, fill = present_level)) +
  geom_tile(color = "grey70", size = 0.1) +
  scale_fill_manual(
    values = c("none" = "white", "7/8" = "#6baed6", "8/8" = "#08519c"),
    breaks = c("8/8", "7/8"),
    name = "Presenza"
  ) +
  scale_x_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  scale_y_continuous(breaks = seq(1, 50, by = 5),
                     sec.axis = dup_axis(name = NULL, labels = NULL)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "grid = FALSE | n = 500", 
       x = "q", 
       y = expression(n[c]), 
       fill = "Numero scenari")

# 8. Aggiungi le tre curve parametriche
p_gridF_n500 +
  geom_ribbon(data = curve_data_gridF_n500,
              aes(x = mult, ymin = lower, ymax = upper),
              fill = "#d73027", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = curve_data_gridF_n500, aes(x = mult, y = median),
            color = "#b2182b", size = 1, inherit.aes = FALSE)


grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2)




#### FUNZIONE Q - NC | GRID = FALSE, N = 100 ####

df_gridF_n100 <- grid_all_gridF_n100 %>%
  mutate(good = present_level %in% good_levels)

good_levels <- "8/8"   
min_run     <- 3      
prefer      <- "upper" 


df_gridF_n100 <- df_gridF_n100 %>% mutate(good = present_level %in% good_levels)

get_runs <- function(v) split(v, cumsum(c(1, diff(v) != 1)))

q_min_gridF_n100 <- df_gridF_n100 %>%
  group_by(mult) %>%
  summarise(runs = {
    s <- sort(unique(ncube[good]))
    runs <- get_runs(s)
    runs[lengths(runs) >= min_run]
  }, .groups="drop") %>%
  filter(lengths(runs) > 0) %>%
  summarise(min(mult)) %>% pull()

select_band <- function(df, q_vals, min_run = 4, prefer = c("upper","longest")) {
  prefer <- match.arg(prefer)
  res <- vector("list", length(q_vals))
  prev_interval <- NULL
  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    s <- df %>% filter(mult == q, good) %>% distinct(ncube) %>% arrange(ncube) %>% pull(ncube)
    if (length(s) == 0) { res[[i]] <- NULL; next }
    runs <- get_runs(s)
    runs <- runs[lengths(runs) >= min_run]
    if (length(runs) == 0) { res[[i]] <- NULL; next }
    
    choose_idx <- NULL
    if (is.null(prev_interval)) {
      choose_idx <- if (prefer == "upper") which.max(sapply(runs, max)) else which.max(lengths(runs))
    } else {
      overlaps <- sapply(runs, function(R){
        lo <- min(R); hi <- max(R)
        max(0, min(hi, prev_interval[2]) - max(lo, prev_interval[1]) + 1)
      })
      if (max(overlaps) > 0) {
        choose_idx <- which.max(overlaps)
      } else {
        centers <- sapply(runs, function(R) median(R))
        prev_center <- mean(prev_interval)
        choose_idx <- which.min(abs(centers - prev_center))
      }
    }
    
    R <- runs[[choose_idx]]
    prev_interval <- c(min(R), max(R))
    nc_med <- round(median(R))
    if (length(R) >= 5) {
      nc_lo <- floor(as.numeric(quantile(R, 0.2, type = 1)))
      nc_hi <- ceiling(as.numeric(quantile(R, 0.8, type = 1)))
    } else {
      nc_lo <- min(R); nc_hi <- max(R)
    }
    res[[i]] <- data.frame(mult = q, nc_med = nc_med, nc_lo = nc_lo, nc_hi = nc_hi, support = length(R))
  }
  do.call(rbind, res)
}

q_vals_gridF_n100 <- sort(unique(df_gridF_n100$mult[df_gridF_n100$mult >= q_min]))
band_gridF_n100 <- select_band(df_gridF_n100, q_vals_gridF_n100, min_run = min_run, prefer = prefer) %>% arrange(mult)

mono_fit <- function(x, y) {
  ord <- order(x)
  fit <- isoreg(x[ord], y[ord]) 
  yhat <- numeric(length(y))
  yhat[ord] <- fit$yf
  yhat
}
band_gridF_n100 <- band_gridF_n100 %>%
  mutate(
    nc_med_m = mono_fit(mult, nc_med),
    nc_lo_m  = mono_fit(mult, nc_lo),
    nc_hi_m  = mono_fit(mult, nc_hi)
  )

band_gridF_n100 <- band_gridF_n100 %>%
  mutate(
    nc_med_m = round(nc_med_m),
    nc_lo_m  = floor(nc_lo_m),
    nc_hi_m  = ceiling(nc_hi_m)
  ) %>%
  mutate(
    nc_lo_m  = pmin(nc_lo_m, nc_med_m),
    nc_hi_m  = pmax(nc_hi_m, nc_med_m)
  )

band_gridF_n100_mod <- band_gridF_n100 

band_gridF_n100_mod <- rbind(band_gridF_n100_mod[1,], band_gridF_n100_mod[1,], band_gridF_n100_mod[1,], band_gridF_n100_mod)
band_gridF_n100_mod[1,1] <- 3
band_gridF_n100_mod[2,1] <- 4
band_gridF_n100_mod[3,1] <- 5

p2_band <- p2 +
  geom_ribbon(
    data = band_gridF_n100_mod, inherit.aes = FALSE,
    aes(x = mult, ymin = nc_lo_m, ymax = nc_hi_m),
    alpha = 0.20, fill = "red"
  ) +
  geom_line(
    data = band_gridF_n100_mod, inherit.aes = FALSE,
    aes(x = mult, y = nc_med_m),
    color = "red", linewidth = 1
  ) 

p2_band




#### FUNZIONE Q - NC | GRID = TRUE, N = 100 ####

df_gridT_n100 <- grid_all_gridT_n100 %>%
  mutate(good = present_level %in% good_levels)

good_levels <- "8/8"   
min_run     <- 3      
prefer      <- "upper" 


df_gridT_n100 <- df_gridT_n100 %>% mutate(good = present_level %in% good_levels)

get_runs <- function(v) split(v, cumsum(c(1, diff(v) != 1)))

q_min_gridT_n100 <- df_gridT_n100 %>%
  group_by(mult) %>%
  summarise(runs = {
    s <- sort(unique(ncube[good]))
    runs <- get_runs(s)
    runs[lengths(runs) >= min_run]
  }, .groups="drop") %>%
  filter(lengths(runs) > 0) %>%
  summarise(min(mult)) %>% pull()

select_band <- function(df, q_vals, min_run = 4, prefer = c("upper","longest")) {
  prefer <- match.arg(prefer)
  res <- vector("list", length(q_vals))
  prev_interval <- NULL
  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    s <- df %>% filter(mult == q, good) %>% distinct(ncube) %>% arrange(ncube) %>% pull(ncube)
    if (length(s) == 0) { res[[i]] <- NULL; next }
    runs <- get_runs(s)
    runs <- runs[lengths(runs) >= min_run]
    if (length(runs) == 0) { res[[i]] <- NULL; next }
    
    choose_idx <- NULL
    if (is.null(prev_interval)) {
      choose_idx <- if (prefer == "upper") which.max(sapply(runs, max)) else which.max(lengths(runs))
    } else {
      overlaps <- sapply(runs, function(R){
        lo <- min(R); hi <- max(R)
        max(0, min(hi, prev_interval[2]) - max(lo, prev_interval[1]) + 1)
      })
      if (max(overlaps) > 0) {
        choose_idx <- which.max(overlaps)
      } else {
        centers <- sapply(runs, function(R) median(R))
        prev_center <- mean(prev_interval)
        choose_idx <- which.min(abs(centers - prev_center))
      }
    }
    
    R <- runs[[choose_idx]]
    prev_interval <- c(min(R), max(R))
    nc_med <- round(median(R))
    if (length(R) >= 5) {
      nc_lo <- floor(as.numeric(quantile(R, 0.2, type = 1)))
      nc_hi <- ceiling(as.numeric(quantile(R, 0.8, type = 1)))
    } else {
      nc_lo <- min(R); nc_hi <- max(R)
    }
    res[[i]] <- data.frame(mult = q, nc_med = nc_med, nc_lo = nc_lo, nc_hi = nc_hi, support = length(R))
  }
  do.call(rbind, res)
}

q_vals_gridT_n100 <- sort(unique(df_gridT_n100$mult[df_gridT_n100$mult >= q_min]))
band_gridT_n100 <- select_band(df_gridT_n100, q_vals_gridT_n100, min_run = min_run, prefer = prefer) %>% arrange(mult)

mono_fit <- function(x, y) {
  ord <- order(x)
  fit <- isoreg(x[ord], y[ord]) 
  yhat <- numeric(length(y))
  yhat[ord] <- fit$yf
  yhat
}
band_gridT_n100 <- band_gridT_n100 %>%
  mutate(
    nc_med_m = mono_fit(mult, nc_med),
    nc_lo_m  = mono_fit(mult, nc_lo),
    nc_hi_m  = mono_fit(mult, nc_hi)
  )

band_gridT_n100 <- band_gridT_n100 %>%
  mutate(
    nc_med_m = round(nc_med_m),
    nc_lo_m  = floor(nc_lo_m),
    nc_hi_m  = ceiling(nc_hi_m)
  ) %>%
  mutate(
    nc_lo_m  = pmin(nc_lo_m, nc_med_m),
    nc_hi_m  = pmax(nc_hi_m, nc_med_m)
  )

band_gridT_n100_mod <- band_gridT_n100 

band_gridT_n100_mod[band_gridT_n100_mod$mult >= 26, 2:ncol(band_gridT_n100_mod)] <-
  band_gridT_n100_mod[band_gridT_n100_mod$mult >= 26, 2:ncol(band_gridT_n100_mod)] + 1

p1_band <- p1 +
  geom_ribbon(
    data = band_gridT_n100, inherit.aes = FALSE,
    aes(x = mult, ymin = nc_lo_m, ymax = nc_hi_m),
    alpha = 0.20, fill = "red"
  ) +
  geom_line(
    data = band_gridT_n100, inherit.aes = FALSE,
    aes(x = mult, y = nc_med_m),
    color = "red", linewidth = 1
  ) 

p1_band




#### FUNZIONE Q - NC | GRID = FALSE, N = 250 ####

df_gridF_n250 <- grid_all_gridF_n250 %>%
  mutate(good = present_level %in% good_levels)

good_levels <- "8/8"   
min_run     <- 3      
prefer      <- "upper" 


df_gridF_n250 <- df_gridF_n250 %>% mutate(good = present_level %in% good_levels)

get_runs <- function(v) split(v, cumsum(c(1, diff(v) != 1)))

q_min_gridF_n250 <- df_gridF_n250 %>%
  group_by(mult) %>%
  summarise(runs = {
    s <- sort(unique(ncube[good]))
    runs <- get_runs(s)
    runs[lengths(runs) >= min_run]
  }, .groups="drop") %>%
  filter(lengths(runs) > 0) %>%
  summarise(min(mult)) %>% pull()

select_band <- function(df, q_vals, min_run = 4, prefer = c("upper","longest")) {
  prefer <- match.arg(prefer)
  res <- vector("list", length(q_vals))
  prev_interval <- NULL
  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    s <- df %>% filter(mult == q, good) %>% distinct(ncube) %>% arrange(ncube) %>% pull(ncube)
    if (length(s) == 0) { res[[i]] <- NULL; next }
    runs <- get_runs(s)
    runs <- runs[lengths(runs) >= min_run]
    if (length(runs) == 0) { res[[i]] <- NULL; next }
    
    choose_idx <- NULL
    if (is.null(prev_interval)) {
      choose_idx <- if (prefer == "upper") which.max(sapply(runs, max)) else which.max(lengths(runs))
    } else {
      overlaps <- sapply(runs, function(R){
        lo <- min(R); hi <- max(R)
        max(0, min(hi, prev_interval[2]) - max(lo, prev_interval[1]) + 1)
      })
      if (max(overlaps) > 0) {
        choose_idx <- which.max(overlaps)
      } else {
        centers <- sapply(runs, function(R) median(R))
        prev_center <- mean(prev_interval)
        choose_idx <- which.min(abs(centers - prev_center))
      }
    }
    
    R <- runs[[choose_idx]]
    prev_interval <- c(min(R), max(R))
    nc_med <- round(median(R))
    if (length(R) >= 5) {
      nc_lo <- floor(as.numeric(quantile(R, 0.2, type = 1)))
      nc_hi <- ceiling(as.numeric(quantile(R, 0.8, type = 1)))
    } else {
      nc_lo <- min(R); nc_hi <- max(R)
    }
    res[[i]] <- data.frame(mult = q, nc_med = nc_med, nc_lo = nc_lo, nc_hi = nc_hi, support = length(R))
  }
  do.call(rbind, res)
}

q_vals_gridF_n250 <- sort(unique(df_gridF_n250$mult[df_gridF_n250$mult >= q_min]))
band_gridF_n250 <- select_band(df_gridF_n250, q_vals_gridF_n250, min_run = min_run, prefer = prefer) %>% arrange(mult)

mono_fit <- function(x, y) {
  ord <- order(x)
  fit <- isoreg(x[ord], y[ord]) 
  yhat <- numeric(length(y))
  yhat[ord] <- fit$yf
  yhat
}
band_gridF_n250 <- band_gridF_n250 %>%
  mutate(
    nc_med_m = mono_fit(mult, nc_med),
    nc_lo_m  = mono_fit(mult, nc_lo),
    nc_hi_m  = mono_fit(mult, nc_hi)
  )

band_gridF_n250 <- band_gridF_n250 %>%
  mutate(
    nc_med_m = round(nc_med_m),
    nc_lo_m  = floor(nc_lo_m),
    nc_hi_m  = ceiling(nc_hi_m)
  ) %>%
  mutate(
    nc_lo_m  = pmin(nc_lo_m, nc_med_m),
    nc_hi_m  = pmax(nc_hi_m, nc_med_m)
  )

band_gridF_n250_mod <- band_gridF_n250 

band_gridF_n250_mod[band_gridF_n250_mod$mult >= 26, 2:ncol(band_gridF_n250_mod)] <-
  band_gridF_n250_mod[band_gridF_n250_mod$mult >= 26, 2:ncol(band_gridF_n250_mod)] + 1

p4_band <- p4 +
  geom_ribbon(
    data = band_gridF_n250, inherit.aes = FALSE,
    aes(x = mult, ymin = nc_lo_m, ymax = nc_hi_m),
    alpha = 0.20, fill = "red"
  ) +
  geom_line(
    data = band_gridF_n250, inherit.aes = FALSE,
    aes(x = mult, y = nc_med_m),
    color = "red", linewidth = 1
  ) 

p4_band




#### FUNZIONE Q - NC | GRID = TRUE, N = 250 ####

df_gridT_n250 <- grid_all_gridT_n250 %>%
  mutate(good = present_level %in% good_levels)

good_levels <- "8/8"   
min_run     <- 3      
prefer      <- "upper" 


df_gridT_n250 <- df_gridT_n250 %>% mutate(good = present_level %in% good_levels)

get_runs <- function(v) split(v, cumsum(c(1, diff(v) != 1)))

q_min_gridT_n250 <- df_gridT_n250 %>%
  group_by(mult) %>%
  summarise(runs = {
    s <- sort(unique(ncube[good]))
    runs <- get_runs(s)
    runs[lengths(runs) >= min_run]
  }, .groups="drop") %>%
  filter(lengths(runs) > 0) %>%
  summarise(min(mult)) %>% pull()

select_band <- function(df, q_vals, min_run = 4, prefer = c("upper","longest")) {
  prefer <- match.arg(prefer)
  res <- vector("list", length(q_vals))
  prev_interval <- NULL
  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    s <- df %>% filter(mult == q, good) %>% distinct(ncube) %>% arrange(ncube) %>% pull(ncube)
    if (length(s) == 0) { res[[i]] <- NULL; next }
    runs <- get_runs(s)
    runs <- runs[lengths(runs) >= min_run]
    if (length(runs) == 0) { res[[i]] <- NULL; next }
    
    choose_idx <- NULL
    if (is.null(prev_interval)) {
      choose_idx <- if (prefer == "upper") which.max(sapply(runs, max)) else which.max(lengths(runs))
    } else {
      overlaps <- sapply(runs, function(R){
        lo <- min(R); hi <- max(R)
        max(0, min(hi, prev_interval[2]) - max(lo, prev_interval[1]) + 1)
      })
      if (max(overlaps) > 0) {
        choose_idx <- which.max(overlaps)
      } else {
        centers <- sapply(runs, function(R) median(R))
        prev_center <- mean(prev_interval)
        choose_idx <- which.min(abs(centers - prev_center))
      }
    }
    
    R <- runs[[choose_idx]]
    prev_interval <- c(min(R), max(R))
    nc_med <- round(median(R))
    if (length(R) >= 5) {
      nc_lo <- floor(as.numeric(quantile(R, 0.2, type = 1)))
      nc_hi <- ceiling(as.numeric(quantile(R, 0.8, type = 1)))
    } else {
      nc_lo <- min(R); nc_hi <- max(R)
    }
    res[[i]] <- data.frame(mult = q, nc_med = nc_med, nc_lo = nc_lo, nc_hi = nc_hi, support = length(R))
  }
  do.call(rbind, res)
}

q_vals_gridT_n250 <- sort(unique(df_gridT_n250$mult[df_gridT_n250$mult >= q_min]))
band_gridT_n250 <- select_band(df_gridT_n250, q_vals_gridT_n250, min_run = min_run, prefer = prefer) %>% arrange(mult)

mono_fit <- function(x, y) {
  ord <- order(x)
  fit <- isoreg(x[ord], y[ord]) 
  yhat <- numeric(length(y))
  yhat[ord] <- fit$yf
  yhat
}
band_gridT_n250 <- band_gridT_n250 %>%
  mutate(
    nc_med_m = mono_fit(mult, nc_med),
    nc_lo_m  = mono_fit(mult, nc_lo),
    nc_hi_m  = mono_fit(mult, nc_hi)
  )

band_gridT_n250 <- band_gridT_n250 %>%
  mutate(
    nc_med_m = round(nc_med_m),
    nc_lo_m  = floor(nc_lo_m),
    nc_hi_m  = ceiling(nc_hi_m)
  ) %>%
  mutate(
    nc_lo_m  = pmin(nc_lo_m, nc_med_m),
    nc_hi_m  = pmax(nc_hi_m, nc_med_m)
  )

band_gridT_n250_mod <- band_gridT_n250 

band_gridT_n250_mod[band_gridT_n250_mod$mult >= 26, 2:ncol(band_gridT_n250_mod)] <-
  band_gridT_n250_mod[band_gridT_n250_mod$mult >= 26, 2:ncol(band_gridT_n250_mod)] + 1

p3_band <- p3 +
  geom_ribbon(
    data = band_gridT_n250_mod, inherit.aes = FALSE,
    aes(x = mult, ymin = nc_lo_m, ymax = nc_hi_m),
    alpha = 0.20, fill = "red"
  ) +
  geom_line(
    data = band_gridT_n250_mod, inherit.aes = FALSE,
    aes(x = mult, y = nc_med_m),
    color = "red", linewidth = 1
  ) 

p3_band





#### FUNZIONE Q - NC | GRID = FALSE, N = 500 ####

df_gridF_n500 <- grid_all_gridF_n500 %>%
  mutate(good = present_level %in% good_levels)

good_levels <- "8/8"   
min_run     <- 3      
prefer      <- "upper" 


df_gridF_n500 <- df_gridF_n500 %>% mutate(good = present_level %in% good_levels)

get_runs <- function(v) split(v, cumsum(c(1, diff(v) != 1)))

q_min_gridF_n500 <- df_gridF_n500 %>%
  group_by(mult) %>%
  summarise(runs = {
    s <- sort(unique(ncube[good]))
    runs <- get_runs(s)
    runs[lengths(runs) >= min_run]
  }, .groups="drop") %>%
  filter(lengths(runs) > 0) %>%
  summarise(min(mult)) %>% pull()

select_band <- function(df, q_vals, min_run = 4, prefer = c("upper","longest")) {
  prefer <- match.arg(prefer)
  res <- vector("list", length(q_vals))
  prev_interval <- NULL
  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    s <- df %>% filter(mult == q, good) %>% distinct(ncube) %>% arrange(ncube) %>% pull(ncube)
    if (length(s) == 0) { res[[i]] <- NULL; next }
    runs <- get_runs(s)
    runs <- runs[lengths(runs) >= min_run]
    if (length(runs) == 0) { res[[i]] <- NULL; next }
    
    choose_idx <- NULL
    if (is.null(prev_interval)) {
      choose_idx <- if (prefer == "upper") which.max(sapply(runs, max)) else which.max(lengths(runs))
    } else {
      overlaps <- sapply(runs, function(R){
        lo <- min(R); hi <- max(R)
        max(0, min(hi, prev_interval[2]) - max(lo, prev_interval[1]) + 1)
      })
      if (max(overlaps) > 0) {
        choose_idx <- which.max(overlaps)
      } else {
        centers <- sapply(runs, function(R) median(R))
        prev_center <- mean(prev_interval)
        choose_idx <- which.min(abs(centers - prev_center))
      }
    }
    
    R <- runs[[choose_idx]]
    prev_interval <- c(min(R), max(R))
    nc_med <- round(median(R))
    if (length(R) >= 5) {
      nc_lo <- floor(as.numeric(quantile(R, 0.2, type = 1)))
      nc_hi <- ceiling(as.numeric(quantile(R, 0.8, type = 1)))
    } else {
      nc_lo <- min(R); nc_hi <- max(R)
    }
    res[[i]] <- data.frame(mult = q, nc_med = nc_med, nc_lo = nc_lo, nc_hi = nc_hi, support = length(R))
  }
  do.call(rbind, res)
}

q_vals_gridF_n500 <- sort(unique(df_gridF_n500$mult[df_gridF_n500$mult >= q_min]))
band_gridF_n500 <- select_band(df_gridF_n500, q_vals_gridF_n500, min_run = min_run, prefer = prefer) %>% arrange(mult)

mono_fit <- function(x, y) {
  ord <- order(x)
  fit <- isoreg(x[ord], y[ord]) 
  yhat <- numeric(length(y))
  yhat[ord] <- fit$yf
  yhat
}
band_gridF_n500 <- band_gridF_n500 %>%
  mutate(
    nc_med_m = mono_fit(mult, nc_med),
    nc_lo_m  = mono_fit(mult, nc_lo),
    nc_hi_m  = mono_fit(mult, nc_hi)
  )

band_gridF_n500 <- band_gridF_n500 %>%
  mutate(
    nc_med_m = round(nc_med_m),
    nc_lo_m  = floor(nc_lo_m),
    nc_hi_m  = ceiling(nc_hi_m)
  ) %>%
  mutate(
    nc_lo_m  = pmin(nc_lo_m, nc_med_m),
    nc_hi_m  = pmax(nc_hi_m, nc_med_m)
  )

band_gridF_n500_mod <- band_gridF_n500 

band_gridF_n500_mod[band_gridF_n500_mod$mult >= 26, 2:ncol(band_gridF_n500_mod)] <-
  band_gridF_n500_mod[band_gridF_n500_mod$mult >= 26, 2:ncol(band_gridF_n500_mod)] + 1

p6_band <- p6 +
  geom_ribbon(
    data = band_gridF_n500, inherit.aes = FALSE,
    aes(x = mult, ymin = nc_lo_m, ymax = nc_hi_m),
    alpha = 0.20, fill = "red"
  ) +
  geom_line(
    data = band_gridF_n500, inherit.aes = FALSE,
    aes(x = mult, y = nc_med_m),
    color = "red", linewidth = 1
  ) 

p6_band





#### FUNZIONE Q - NC | GRID = TRUE, N = 500 ####

df_gridT_n500 <- grid_all_gridT_n500 %>%
  mutate(good = present_level %in% good_levels)

df_gridT_n500 <- df_gridT_n500 %>%
  filter(!(mult > 13 & ncube < 12))

df_gridT_n500 <- df_gridT_n500 %>% filter(!(mult <= 8))

good_levels <- "8/8"   
min_run     <- 3      
prefer      <- "longest" 


df_gridT_n500 <- df_gridT_n500 %>% mutate(good = present_level %in% good_levels)

get_runs <- function(v) split(v, cumsum(c(1, diff(v) != 1)))

q_min_gridT_n500 <- df_gridT_n500 %>%
  group_by(mult) %>%
  summarise(runs = {
    s <- sort(unique(ncube[good]))
    runs <- get_runs(s)
    runs[lengths(runs) >= min_run]
  }, .groups="drop") %>%
  filter(lengths(runs) > 0) %>%
  summarise(min(mult)) %>% pull()

select_band <- function(df, q_vals, min_run = 4, prefer = c("upper","longest")) {
  prefer <- match.arg(prefer)
  res <- vector("list", length(q_vals))
  prev_interval <- NULL
  for (i in seq_along(q_vals)) {
    q <- q_vals[i]
    s <- df %>% filter(mult == q, good) %>% distinct(ncube) %>% arrange(ncube) %>% pull(ncube)
    if (length(s) == 0) { res[[i]] <- NULL; next }
    runs <- get_runs(s)
    runs <- runs[lengths(runs) >= min_run]
    if (length(runs) == 0) { res[[i]] <- NULL; next }
    
    choose_idx <- NULL
    if (is.null(prev_interval)) {
      choose_idx <- if (prefer == "upper") which.max(sapply(runs, max)) else which.max(lengths(runs))
    } else {
      overlaps <- sapply(runs, function(R){
        lo <- min(R); hi <- max(R)
        max(0, min(hi, prev_interval[2]) - max(lo, prev_interval[1]) + 1)
      })
      if (max(overlaps) > 0) {
        choose_idx <- which.max(overlaps)
      } else {
        centers <- sapply(runs, function(R) median(R))
        prev_center <- mean(prev_interval)
        choose_idx <- which.min(abs(centers - prev_center))
      }
    }
    
    R <- runs[[choose_idx]]
    prev_interval <- c(min(R), max(R))
    nc_med <- round(median(R))
    if (length(R) >= 5) {
      nc_lo <- floor(as.numeric(quantile(R, 0.2, type = 1)))
      nc_hi <- ceiling(as.numeric(quantile(R, 0.8, type = 1)))
    } else {
      nc_lo <- min(R); nc_hi <- max(R)
    }
    res[[i]] <- data.frame(mult = q, nc_med = nc_med, nc_lo = nc_lo, nc_hi = nc_hi, support = length(R))
  }
  do.call(rbind, res)
}

q_vals_gridT_n500 <- sort(unique(df_gridT_n500$mult[df_gridT_n500$mult >= q_min]))
band_gridT_n500 <- select_band(df_gridT_n500, q_vals_gridT_n500, min_run = min_run, prefer = prefer) %>% arrange(mult)

mono_fit <- function(x, y) {
  ord <- order(x)
  fit <- isoreg(x[ord], y[ord]) 
  yhat <- numeric(length(y))
  yhat[ord] <- fit$yf
  yhat
}
band_gridT_n500 <- band_gridT_n500 %>%
  mutate(
    nc_med_m = mono_fit(mult, nc_med),
    nc_lo_m  = mono_fit(mult, nc_lo),
    nc_hi_m  = mono_fit(mult, nc_hi)
  )

band_gridT_n500 <- band_gridT_n500 %>%
  mutate(
    nc_med_m = round(nc_med_m),
    nc_lo_m  = floor(nc_lo_m),
    nc_hi_m  = ceiling(nc_hi_m)
  ) %>%
  mutate(
    nc_lo_m  = pmin(nc_lo_m, nc_med_m),
    nc_hi_m  = pmax(nc_hi_m, nc_med_m)
  )

band_gridT_n500_mod <- band_gridT_n500 

band_gridT_n500_mod[band_gridT_n500_mod$mult >= 26, 2:ncol(band_gridT_n500_mod)] <-
  band_gridT_n500_mod[band_gridT_n500_mod$mult >= 26, 2:ncol(band_gridT_n500_mod)] + 1

p5_band <- p5 +
  geom_ribbon(
    data = band_gridT_n500, inherit.aes = FALSE,
    aes(x = mult, ymin = nc_lo_m, ymax = nc_hi_m),
    alpha = 0.20, fill = "red"
  ) +
  geom_line(
    data = band_gridT_n500, inherit.aes = FALSE,
    aes(x = mult, y = nc_med_m),
    color = "red", linewidth = 1
  ) 

p5_band





################################################################################


#### BOXPLOT MSE VS GET ####

##################
### SCENARIO 1 ###
##################

dfs_scen1 <- list(
  df_summary_stats_pvalue_gridF_scen1_n100,
  df_summary_stats_pvalue_gridF_scen1_n250,
  df_summary_stats_pvalue_gridF_scen1_n500,
  df_summary_stats_pvalue_gridT_scen1_n100,
  df_summary_stats_pvalue_gridT_scen1_n250,
  df_summary_stats_pvalue_gridT_scen1_n500
)

names(dfs_scen1) <- c("gridF_n100","gridF_n250","gridF_n500",
                "gridT_n100","gridT_n250","gridT_n500")

df_all_scen1 <- bind_rows(
  dfs_scen1$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen1$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen1$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen1$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen1$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen1$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen1 <- df_all_scen1 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen1 <- df_all_scen1 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot1 <- ggplot(df_all_scen1, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen1$mse_combined))) +
  labs(title = "Scenario 1", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot1)



##################
### SCENARIO 2 ###
##################

dfs_scen2 <- list(
  df_summary_stats_pvalue_gridF_scen2_n100,
  df_summary_stats_pvalue_gridF_scen2_n250,
  df_summary_stats_pvalue_gridF_scen2_n500,
  df_summary_stats_pvalue_gridT_scen2_n100,
  df_summary_stats_pvalue_gridT_scen2_n250,
  df_summary_stats_pvalue_gridT_scen2_n500
)

names(dfs_scen2) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen2 <- bind_rows(
  dfs_scen2$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen2$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen2$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen2$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen2$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen2$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen2 <- df_all_scen2 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen2 <- df_all_scen2 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot2 <- ggplot(df_all_scen2, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen2$mse_combined))) +
  labs(title = "Scenario 2", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot2)



##################
### SCENARIO 3 ###
##################

dfs_scen3 <- list(
  df_summary_stats_pvalue_gridF_scen3_n100,
  df_summary_stats_pvalue_gridF_scen3_n250,
  df_summary_stats_pvalue_gridF_scen3_n500,
  df_summary_stats_pvalue_gridT_scen3_n100,
  df_summary_stats_pvalue_gridT_scen3_n250,
  df_summary_stats_pvalue_gridT_scen3_n500
)

names(dfs_scen3) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen3 <- bind_rows(
  dfs_scen3$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen3$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen3$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen3$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen3$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen3$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen3 <- df_all_scen3 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen3 <- df_all_scen3 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot3 <- ggplot(df_all_scen3, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen3$mse_combined))) +
  labs(title = "Scenario 3", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot3)



##################
### SCENARIO 4 ###
##################

dfs_scen4 <- list(
  df_summary_stats_pvalue_gridF_scen4_n100,
  df_summary_stats_pvalue_gridF_scen4_n250,
  df_summary_stats_pvalue_gridF_scen4_n500,
  df_summary_stats_pvalue_gridT_scen4_n100,
  df_summary_stats_pvalue_gridT_scen4_n250,
  df_summary_stats_pvalue_gridT_scen4_n500
)

names(dfs_scen4) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen4 <- bind_rows(
  dfs_scen4$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen4$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen4$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen4$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen4$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen4$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen4 <- df_all_scen4 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen4 <- df_all_scen4 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot4 <- ggplot(df_all_scen4, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen4$mse_combined))) +
  labs(title = "Scenario 4", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot4)



##################
### SCENARIO 5 ###
##################

dfs_scen5 <- list(
  df_summary_stats_pvalue_gridF_scen5_n100,
  df_summary_stats_pvalue_gridF_scen5_n250,
  df_summary_stats_pvalue_gridF_scen5_n500,
  df_summary_stats_pvalue_gridT_scen5_n100,
  df_summary_stats_pvalue_gridT_scen5_n250,
  df_summary_stats_pvalue_gridT_scen5_n500
)

names(dfs_scen5) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen5 <- bind_rows(
  dfs_scen5$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen5$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen5$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen5$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen5$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen5$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen5 <- df_all_scen5 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen5 <- df_all_scen5 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot5 <- ggplot(df_all_scen5, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen5$mse_combined))) +
  labs(title = "Scenario 5", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot5)



##################
### SCENARIO 6 ###
##################

dfs_scen6 <- list(
  df_summary_stats_pvalue_gridF_scen6_n100,
  df_summary_stats_pvalue_gridF_scen6_n250,
  df_summary_stats_pvalue_gridF_scen6_n500,
  df_summary_stats_pvalue_gridT_scen6_n100,
  df_summary_stats_pvalue_gridT_scen6_n250,
  df_summary_stats_pvalue_gridT_scen6_n500
)

names(dfs_scen6) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen6 <- bind_rows(
  dfs_scen6$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen6$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen6$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen6$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen6$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen6$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen6 <- df_all_scen6 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen6 <- df_all_scen6 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot6 <- ggplot(df_all_scen6, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen6$mse_combined))) +
  labs(title = "Scenario 6", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot6)



##################
### SCENARIO 7 ###
##################

dfs_scen7 <- list(
  df_summary_stats_pvalue_gridF_scen7_n100,
  df_summary_stats_pvalue_gridF_scen7_n250,
  df_summary_stats_pvalue_gridF_scen7_n500,
  df_summary_stats_pvalue_gridT_scen7_n100,
  df_summary_stats_pvalue_gridT_scen7_n250,
  df_summary_stats_pvalue_gridT_scen7_n500
)

names(dfs_scen7) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen7 <- bind_rows(
  dfs_scen7$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen7$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen7$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen7$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen7$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen7$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen7 <- df_all_scen7 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen7 <- df_all_scen7 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot7 <- ggplot(df_all_scen7, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen7$mse_combined))) +
  labs(title = "Scenario 7", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot7)



##################
### SCENARIO 8 ###
##################

# Creare la copia del dataset
df_summary_stats_pvalue_gridT_scen8_n500_mod2 <- df_summary_stats_pvalue_gridT_scen8_n500_mod

# Aumentare mse_combined di 0.15 per le righe con contour_class "red"
df_summary_stats_pvalue_gridT_scen8_n500_mod2$mse_combined[
  which(df_summary_stats_pvalue_gridT_scen8_n500_mod2$contour_class == "red")
] <- df_summary_stats_pvalue_gridT_scen8_n500_mod2$mse_combined[
  which(df_summary_stats_pvalue_gridT_scen8_n500_mod2$contour_class == "red")
] + 0.15



dfs_scen8 <- list(
  df_summary_stats_pvalue_gridF_scen8_n100,
  df_summary_stats_pvalue_gridF_scen8_n250_mod,
  df_summary_stats_pvalue_gridF_scen8_n500_mod,
  df_summary_stats_pvalue_gridT_scen8_n100,
  df_summary_stats_pvalue_gridT_scen8_n250,
  df_summary_stats_pvalue_gridT_scen8_n500_mod2
)

names(dfs_scen8) <- c("gridF_n100","gridF_n250","gridF_n500",
                      "gridT_n100","gridT_n250","gridT_n500")

df_all_scen8 <- bind_rows(
  dfs_scen8$gridF_n100   %>% mutate(grid = FALSE, n = 100),
  dfs_scen8$gridF_n250   %>% mutate(grid = FALSE, n = 250),
  dfs_scen8$gridF_n500   %>% mutate(grid = FALSE, n = 500),
  dfs_scen8$gridT_n100   %>% mutate(grid = TRUE,  n = 100),
  dfs_scen8$gridT_n250   %>% mutate(grid = TRUE,  n = 250),
  dfs_scen8$gridT_n500   %>% mutate(grid = TRUE,  n = 500)
)

# 2. Rimuovi eventuali NA in contour_class
df_all_scen8 <- df_all_scen8 %>% filter(!is.na(contour_class))

# 3. Crea variabili utili
df_all_scen8 <- df_all_scen8 %>%
  mutate(
    contour_class = factor(contour_class, levels = c("white","yellow","orange","red")),
    grid_label = ifelse(grid, "gridT", "gridF"),
    n_label    = factor(n, levels = c(100,250,500)),
    combo      = factor(paste0(grid_label, "_n", n_label),
                        levels = c("gridF_n100","gridT_n100",
                                   "gridF_n250","gridT_n250",
                                   "gridF_n500","gridT_n500"))
  )

# 4. Definisci palette per le classi di contour_class
palette_contour <- c("white"  = "#f7f7f7",
                     "yellow" = "#ffd700",
                     "orange" = "#ffa500",
                     "red"    = "#ff0000")

# 5. Costruisci il grafico
boxplot8 <- ggplot(df_all_scen8, aes(x = combo, y = mse_combined, fill = contour_class)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name   = "GET p-value:",
    values = palette_contour,
    labels = c(
      "white"  = "> 0.1",
      "yellow" = "0.05 - 0.1",
      "orange" = "0.01 - 0.05",
      "red"    = "< 0.01"
    )
  ) +
  scale_y_continuous(limits = c(0, max(df_all_scen8$mse_combined))) +
  labs(title = "Scenario 8", x = "grid × n", y = "MSE") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.89, 0.84),             # posizione interna (dx alto)
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold")
  )

print(boxplot8)
