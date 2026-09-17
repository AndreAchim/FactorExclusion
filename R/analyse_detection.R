# analyse_detection.R
# Ajuste une hierarchie de quatre modeles GLM quasibinomial emboites et compare
# par tests F pour documenter l'effet des saturations des FIVs (sat) sur la
# puissance de detection des correlations residuelles.
#
# Prerequisite : df_det doit exister dans l'environnement (produit par
#   build_detection_data.R).
#
# Usage :
#   res <- analyse_detection()          # utilise df_det global
#   res <- analyse_detection(df_det)    # explicite
#
# Valeur de retour : liste invisible avec $m0..$m3 et $anova

analyse_detection <- function(dat = df_det) {

  # m0 : reference — att et N seulement, sans sat
  m0 <- glm(cbind(n_det, n_tot - n_det) ~ log(att) * log(N_val),
             data = dat, family = quasibinomial)

  # m1 : effet principal de sat (meme pente pour att et N dans tous les groupes)
  m1 <- glm(cbind(n_det, n_tot - n_det) ~ log(att) * log(N_val) + sat,
             data = dat, family = quasibinomial)

  # m2 : interactions 2-voies sat x log(att) et sat x log(N_val)
  #       (equivalent a m3q de build_detection_data.R)
  m2 <- glm(cbind(n_det, n_tot - n_det) ~
              log(att) * sat + log(N_val) * sat + log(att):log(N_val),
             data = dat, family = quasibinomial)

  # m3 : interaction triple — sat modifie aussi l'effet conjoint att x N
  m3 <- glm(cbind(n_det, n_tot - n_det) ~ log(att) * log(N_val) * sat,
             data = dat, family = quasibinomial)

  # Tests F quasi emboites
  av <- anova(m0, m1, m2, m3, test = "F")

  cat("\n=== Effet de sat sur la puissance de detection ===\n")
  cat("Modeles quasibinomiaux emboites, tests F\n\n")
  cat("  m0 : ~ log(att) * log(N_val)               [reference, sans sat]\n")
  cat("  m1 : m0 + sat                               [effet principal]\n")
  cat("  m2 : ~ log(att)*sat + log(N_val)*sat + ...  [interactions 2-voies]\n")
  cat("  m3 : ~ log(att) * log(N_val) * sat          [interaction triple]\n\n")
  print(av)

  # Dispersion estimee (basee sur m3, modele le plus sature)
  phi <- summary(m3)$dispersion
  cat(sprintf("\nDispersion estimee (m3) : %.1f\n", phi))

  invisible(list(m0 = m0, m1 = m1, m2 = m2, m3 = m3, anova = av))
}
