partitionne_R <- function(AS, ote, seuil = .05) {
  # Estime les saturations et les corrélations résiduelles pour toutes
  # les variables hors ancres, à partir de la structure SCA initialisée.
  R <- AS$R
  N <- AS$N
  t_crit <- qt(1 - seuil / 2, df = N - 2)
  r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  if (abs(R[ote[1], ote[2]]) < r_crit) {
    return(NULL)
  }
  nv    <- nrow(R)
  sat   <- sat_sur_facteur(R, N, ote)
  sa    <- sat$details$x_est
  Heywood <- any(abs(sa) > 1)
  if (Heywood) warning("Condition de Heywood détectée.")
  SA      <- sa %*% t(sa)
  p2      <- R[-ote, -ote] - SA
  diag(p2) <- NA
  p2[upper.tri(p2)] <- prob_corr_residuelles(R, ote, N)
  return(list(corr_resid = round(p2, 4), details = sat$details, Heywood = Heywood))
}
