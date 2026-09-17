partitionne_R <- function(AS,ote,seuil=.05){
  # on aura commencé, en produisant AS, par tester l'annulation des variables
  R <- AS$R
  N <- AS$N
  # Valeur critique de |r| pour ce N
  t_crit <- qt(1 - seuil / 2, df = N - 2)
  r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  # valider que les covariables sont corrélées
  if (abs(R[ote[1],ote[2]]) < r_crit){
    # cat("Correlation of",R[ote[1],ote[2]],"is below critical value of", r_crit,"\n")
    return(NULL)
  }
  nv <- nrow(R)
  sat <- sat_sur_facteur(R, N, ote)
  # if (any(sat$sat>1)) browser()
  sa <- sat$details$x_est
  Heywood <- any(abs(sa)>1)
  if (Heywood){
    # browser()
    warning("Heywood condition detected, suggesting unreliable data!\n")
  }
  SA <- sa %*% t(sa)
  produit <- SA[upper.tri(SA)]
  correl <- R[-ote,-ote]
  p2 <- correl- SA
  correl <- correl[upper.tri(correl)]
  # browser()
  diag(p2) <- NA
  p2[upper.tri(p2)] <- prob_corr_residuelles(R, ote, N)
  return(list(corr_resid=round(p2,4),details=sat$details,Heywood=Heywood))
}