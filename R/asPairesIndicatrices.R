asPairesIndicatrices <- function(AS) {
  # Identifier les variables indicatrices, s'il y en a pour la branche AS$branche
  AS$Cpaires <- combn(AS$pertinent, 2) # tous les sous-ensembles de 2 parmi les variables pertinentes
  nc <- ncol(AS$Cpaires)
  AS$Crit <- rep(0, nc)
  AS$Prob <- rep(0, nc)
  AS$Corr <- matrix(0, AS$nv, nc)
  AS$Ppaires <- rep(0, nc)
  AS$doublet <- NULL
  for (k in 1:nc) {
    AS <- asAnnulePaire(AS, k)  # optimiser les variables de la rangée k de AS$Cpaires
  }
  AS$z2 <- (AS$Corr)^2 * (AS$N - 1)
  return(AS)
}
