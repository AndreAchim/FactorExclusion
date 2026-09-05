sc1 <- function(x) {
  # Normalise x à une somme de carrés de 1, sauf si tous des 0
  flg <- any(abs(x) > 1e-10)
  if (flg)
    x <- x / sqrt(sum(x^2))
  return(x)
}

SC1 <- function(M) {
  # Comme sc1 mais pour les colonnes d'une matrice M
  MM <- M
  for (col in 1:ncol(M))
    MM[, col] <- sc1(M[, col])
  return(MM)
}
