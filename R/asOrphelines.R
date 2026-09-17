asOrphelines <- function(AS) {
  # utiliser AS$R pour identifier ses colonnes dont la plus grande valeur
  # absolue (excluant la diagonale) n'est pas significative (p > .01)
  R <- AS$R - diag(AS$nv)
  mxR <- apply(abs(R), 2, max)
  z2 <- mxR^2 * (AS$N - 1)
  AS$pOrpheline <- 1 - pchisq(z2, 1)^(AS$nv - 1)
  # f <- which((AS$pOrpheline > .01) | (mxR < .1))
  f <- which(AS$pOrpheline > AS$seuils[1])  # POC: Why seuils[1]?
  AS$orphelines <- f
  if (length(f) > 0) {
    AS$colOrphelines <- as.matrix(AS$GS[, f])
    AS$GS[, f] <- 0 # pour rendre toutes les corrélations nulles sans mélanger les rangs des variables
  }
  AS$pertinent <- setdiff(1:AS$nv, AS$orphelines)
  return(AS)
}
