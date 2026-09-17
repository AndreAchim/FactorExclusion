# fig_biais.R — Figure de biais des corrélations résiduelles simulées
#
# Prérequis : les objets suivants doivent exister dans l'environnement :
#   efg       — liste de 16 éléments produite par FE_tests()
#   CRatt_all — efg[[16]]$CRatt  (66 valeurs attendues de population)
#   zro       — efg[[16]]$zro    (indices des 20 paires à valeur attendue nulle)
#
# Usage :
#   CRatt_all <- efg[[16]]$CRatt
#   zro       <- efg[[16]]$zro
#   source("R/fig_biais.R")   # produit Figure1.pdf et Figure1.png

library(ggplot2)

# Lire CRatt_all et zro depuis efg pour garantir la cohérence des dimensions
CRatt_all <- efg[[16]]$CRatt
zro       <- efg[[16]]$zro

# --- Ordre global des 46 paires non-nulles (commun à tous les panneaux) ---
# Tri primaire  : |valeur attendue| croissant
# Tri secondaire : dans chaque groupe ex æquo, moyenne observée |moy|
#                  sur les 15 conditions croissante
nonzro_all <- setdiff(seq_along(CRatt_all), zro)   # 46 indices dans CRatt_all
att_signed  <- CRatt_all[nonzro_all]                # valeurs signées
att_abs_all <- abs(att_signed)
s_all       <- sign(att_signed)

grand_moy <- Reduce("+", lapply(seq_len(15), function(i) {
  d4 <- efg[[i]]$stt$CR[4, nonzro_all]
  tanh(atanh(d4) + atanh(att_signed)) * s_all
})) / 15

# global_o2[k] = indice dans nonzro_all de la k-ième paire en ordre d'affichage
global_o2 <- order(att_abs_all, grand_moy)

# --- Extraction des données pour une condition i (1 à 15) ---
extraire_cond <- function(i) {
  cond <- efg[[i]]
  cr   <- cond$stt$CR          # matrice 6 × 66 produite par sim_stats()

  # Indices dans CRatt_all, dans l'ordre global
  idx  <- nonzro_all[global_o2]
  att0 <- CRatt_all[idx]       # valeurs signées
  s    <- sign(att0)

  # Lignes de cr : 4 = tanh(déviation moy), 5 = tanh(IC bas), 6 = tanh(IC haut)
  d4 <- cr[4, idx]; d5 <- cr[5, idx]; d6 <- cr[6, idx]
  df_t <- cr[2, idx]           # degrés de liberté (= nOK - 1)

  moy_atanh  <- atanh(d4) + atanh(att0)
  se_atanh   <- (atanh(d6) - atanh(d5)) / (2 * qt(0.975, df_t))
  demi_ci999 <- qt(0.9995, df_t) * se_atanh
  sd_atanh   <- se_atanh * sqrt(df_t + 1)

  lower_ci <- tanh(moy_atanh - demi_ci999)
  upper_ci <- tanh(moy_atanh + demi_ci999)
  lower_sd <- tanh(moy_atanh - sd_atanh)
  upper_sd <- tanh(moy_atanh + sd_atanh)

  # Inversion de signe pour les paires à valeur attendue négative
  moy_abs    <- tanh(moy_atanh) * s
  lower_ci_p <- ifelse(s > 0, lower_ci, -upper_ci)
  upper_ci_p <- ifelse(s > 0, upper_ci, -lower_ci)
  lower_sd_p <- ifelse(s > 0, lower_sd, -upper_sd)
  upper_sd_p <- ifelse(s > 0, upper_sd, -lower_sd)

  data.frame(
    x        = seq_along(idx),
    att      = abs(att0),
    moy      = moy_abs,
    lower_ci = lower_ci_p, upper_ci = upper_ci_p,
    lower_sd = lower_sd_p, upper_sd = upper_sd_p,
    sat_lab  = paste0("(", cond$sat[1], ", ", cond$sat[2], ")"),
    N_lab    = paste0("N = ", abs(cond$ns)),
    sat1     = cond$sat[1],
    N_val    = abs(cond$ns)
  )
}

# --- Construction du data.frame combiné (15 conditions × 46 paires) ---
df_all <- do.call(rbind, lapply(seq_len(15), extraire_cond))

df_all$sat_lab <- factor(df_all$sat_lab,
  levels = unique(df_all$sat_lab[order(df_all$sat1)]))
df_all$N_lab <- factor(df_all$N_lab,
  levels = paste0("N = ", c(125, 500, 2000)))

# --- Figure 15 panneaux (5 rangées × 3 colonnes) ---
fig15 <- ggplot(df_all, aes(x = x)) +
  geom_ribbon(aes(ymin = lower_sd, ymax = upper_sd),
              fill = "grey60", color = NA) +
  geom_line(aes(y = att, linewidth = "Expected")) +
  geom_segment(
    data      = \(d) d[d$att < d$lower_ci | d$att > d$upper_ci, ],
    aes(x = x, xend = x, y = lower_ci, yend = upper_ci),
    linewidth = 0.25, color = "grey85"
  ) +
  geom_line(aes(y = moy, linewidth = "Mean observed"), color = "white") +
  scale_linewidth_manual(
    name   = NULL,
    values = c("Expected" = 0.5, "Mean observed" = 0.25)
  ) +
  facet_grid(sat_lab ~ N_lab) +
  labs(x = "Variable pairs (ordered by |expected value|, then mean observed within ties)",
       y = "Residual correlation") +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.key       = element_rect(fill = "grey60"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text.x     = element_text(colour = "black"),
        strip.text.y     = element_text(colour = "black", angle = -90))

ggsave("Figure1.pdf", fig15, width = 17, height = 22, units = "cm")
ggsave("Figure1.png", fig15, width = 17, height = 22, units = "cm", dpi = 300)
message("Figure1.pdf et Figure1.png enregistrés.")
