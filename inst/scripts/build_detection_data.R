# build_detection_data.R
# Construit df_det et ajuste les modeles GLM quasibinomial (m0q, m3q).
#
# Prerequisite : efg doit exister dans l'environnement.
#   efg <- FE_tests()   # simulations completes (~15 min)
#   # ou recharger depuis un fichier sauvegarde :
#   efg <- readRDS("efg.rds")
#
# Usage :
#   source("R/build_detection_data.R")
#   # puis sourcer R/fig_detection.R pour produire Figure2.pdf
#
# Structure de efg[[i]] (i = 1..15) :
#   $sat, $ns, $stt, $sigCR (logique nOK x 66), $OK_idx
#   $stt$sgCR : matrice 2 x 66 — ligne 1 = taux de detection a alpha=.05
#   $stt$nOK  : nombre de simulations valides (sans cas Heywood)
# efg[[16]] : $CRatt (66 valeurs attendues), $zro (indices des paires nulles)

# Lire depuis efg pour garantir la cohérence des dimensions
CRatt_all <- efg[[16]]$CRatt
zro       <- efg[[16]]$zro
nonzro    <- setdiff(seq_along(CRatt_all), zro)

# --- Construction de df_det (15 conditions x 46 paires = 690 lignes) ---
df_det <- do.call(rbind, lapply(seq_len(15), function(i) {
  cond <- efg[[i]]
  nOK  <- cond$stt$nOK
  # Utiliser sigCR si disponible (comptage exact), sinon reconstituer via stt
  if (!is.null(cond$sigCR)) {
    n_det <- colSums(cond$sigCR[, nonzro])
  } else {
    n_det <- round(cond$stt$sgCR[1, nonzro] * nOK)
  }
  data.frame(
    n_det = n_det,
    n_tot = nOK,
    att   = abs(CRatt_all[nonzro]),
    N_val = abs(cond$ns),
    sat   = paste0("(", cond$sat[1], ",", cond$sat[2], ")")
  )
}))
df_det$sat <- factor(df_det$sat)

# --- Modeles GLM quasibinomial ---
m0q <- glm(cbind(n_det, n_tot - n_det) ~ log(att) * log(N_val),
           data = df_det, family = quasibinomial)

m3q <- glm(cbind(n_det, n_tot - n_det) ~
             log(att) * sat + log(N_val) * sat + log(att):log(N_val),
           data = df_det, family = quasibinomial)

message("df_det (", nrow(df_det), " lignes) et m0q/m3q crees.")
