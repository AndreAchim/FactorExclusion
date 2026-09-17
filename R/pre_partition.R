pre_partition <- function(R,N,ote,seuil=.05){
  # Valeur critique de |r| — calculée une seule fois pour ce N
  t_crit <- qt(1 - seuil / 2, df = N - 2)
  r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  # valider que les covariables sont corrélées
  if (abs(R[ote[1],ote[2]]) < r_crit){
    warning("Correlation of",R[ote[1],ote[2]],"is below critical value of", r_crit,"\n")
    return(NULL)
  }
  # obtenir toutes les annualtions par paires de variables
  if (is.null(colnames(R)))
    rownames(R) <- colnames(R) <- letters[1:ncol(R)]
  AS <- init_SCA(R,N)
  AS$ote <- ote
  return(AS)   
}