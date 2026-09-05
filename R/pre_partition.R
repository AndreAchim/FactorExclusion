pre_partition <- function(R, N, ote, seuil = .05) {
  # Vérifie que les ancres sont corrélées et initialise la structure SCA.
  t_crit <- qt(1 - seuil / 2, df = N - 2)
  r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  if (abs(R[ote[1], ote[2]]) < r_crit) {
    warning("Corrélation entre ancres inférieure à la valeur critique.")
    return(NULL)
  }
  if (is.null(colnames(R)))
    rownames(R) <- colnames(R) <- letters[1:ncol(R)]
  AS      <- init_SCA(R, N)
  AS$ote  <- ote
  return(AS)
}
