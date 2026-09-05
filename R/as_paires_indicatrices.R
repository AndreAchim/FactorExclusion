as_paires_indicatrices <- function(AS) {
  # Identifier les variables indicatrices en vérifiant l'annulation de
  # toutes les paires de variables non-orphelines.
  AS$Cpaires  <- combn(AS$pertinent, 2)
  nc          <- ncol(AS$Cpaires)
  AS$Crit     <- rep(0, nc)
  AS$Prob     <- rep(0, nc)
  AS$Ppaires  <- rep(0, nc)
  AS$satPaires <- matrix(NA, nrow = 2, ncol = nc)
  AS$doublet  <- NULL
  ASm         <- matrix(0, nrow = AS$nv, ncol = AS$nv)
  rownames(ASm) <- colnames(ASm) <- colnames(AS$R)
  for (k in 1:nc) {
    AS <- as_annule_paire(AS, k)
    ASm[AS$Cpaires[1, k], AS$Cpaires[2, k]] <- AS$Prob[k]
    ASm[AS$Cpaires[2, k], AS$Cpaires[1, k]] <- AS$Ppaires[k]
  }
  diag(ASm) <- NA
  AS$mPaires <- ASm
  return(AS)
}
