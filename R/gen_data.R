gen_data <- function(F, N, R = NULL, with_var = FALSE, ON = FALSE) {
  # Génère N observations à partir d'une matrice de saturations F.
  #
  # Arguments :
  #   F        : matrice de saturations (nv × nf), ou vecteur si nf = 1
  #   N        : taille d'échantillon
  #   R        : matrice de corrélation des facteurs (optionnel)
  #   with_var : si TRUE, la dernière colonne de F contient les variances d'unicité
  #   ON       : si TRUE, génère des données avec corrélations empiriques exactes
  if (!is.matrix(F)) {
    nf   <- 1
    nv   <- length(F)
    unic <- 1 - F * F
    f    <- as.matrix(F)
  } else {
    nv <- nrow(F)
    nf <- ncol(F)
    if (with_var) {
      unic <- F[, nf]
      F    <- F[, -nf]
      nf   <- nf - 1
    } else {
      unic <- 1 - rowSums(F * F)
    }
  }

  if (!is.null(R)) {
    if (isSymmetric(R))
      R <- t(chol(R))
    F <- F %*% R
  }

  src_cov <- diag(nf + nv)
  G       <- cbind(F, diag(sqrt(unic)))
  srce    <- MASS::mvrnorm(N, rep(0, nf + nv), src_cov, empirical = ON)
  dat     <- srce %*% t(G)
  colnames(dat) <- letters[1:nv]
  return(list(dt = dat, srce = srce))
}
