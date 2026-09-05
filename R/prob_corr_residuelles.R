# Covariance asymptotique entre deux corrélations de Pearson r_ij et r_kl
# estimées dans un même échantillon de taille N (formule de Pearson-Filon).
.cov_r <- function(R, i, j, k, l, N) {
  rij <- R[i, j]; rkl <- R[k, l]
  rik <- R[i, k]; ril <- R[i, l]
  rjk <- R[j, k]; rjl <- R[j, l]
  (
    rik * rjl + ril * rjk
    - rij * (rik * ril + rjk * rjl)
    - rkl * (rik * rjk + ril * rjl)
    + (rij * rkl / 2) * (rik^2 + ril^2 + rjk^2 + rjl^2)
  ) / (N - 1)
}

prob_corr_residuelles <- function(R, AB, N) {
  # Test bilatéral de H0 : r(X,Y) = xy  (pas de corrélation résiduelle)
  # par la méthode delta sur atanh(rXY) - atanh(xy), où
  #   xy = saturation(X) × saturation(Y) estimées via les ancres A et B.
  #
  # Arguments :
  #   R  : matrice de corrélations complète (nv × nv)
  #   AB : indices des deux ancres (vecteur de longueur 2)
  #   N  : taille d'échantillon
  #
  # Retourne un vecteur de p-valeurs dans l'ordre du triangle supérieur
  # de R[-AB, -AB] (compatible avec upper.tri, ordre colonne-majeure).

  ia <- AB[1]; ib <- AB[2]
  rAB <- R[ia, ib]

  autres    <- setdiff(seq_len(nrow(R)), AB)
  n_a       <- length(autres)
  idx_pairs <- which(upper.tri(matrix(0, n_a, n_a)), arr.ind = TRUE)
  np        <- nrow(idx_pairs)
  pvals     <- numeric(np)

  for (s in seq_len(np)) {
    ix <- autres[idx_pairs[s, 1]]
    iy <- autres[idx_pairs[s, 2]]

    rAX <- R[ia, ix]; rBX <- R[ib, ix]
    rAY <- R[ia, iy]; rBY <- R[ib, iy]
    rXY <- R[ix, iy]

    x2 <- rAX * rBX / rAB
    y2 <- rAY * rBY / rAB
    x2 <- max(x2, 1e-6)
    y2 <- max(y2, 1e-6)

    xy <- sqrt(x2 * y2) * sign(rAX) * sign(rAY)
    if (abs(xy) >= 1 || abs(rXY) >= 1) { pvals[s] <- NA; next }

    dz <- atanh(rXY) - atanh(xy)

    # Gradient de (atanh(rXY) - atanh(xy)) par rapport aux 6 corrélations
    # Ordre : rXY, rAX, rBX, rAY, rBY, rAB
    c2xy  <- 1 - xy^2
    c2rXY <- 1 - rXY^2
    g <- c(
       1 / c2rXY,
      -xy / (2 * rAX * c2xy),
      -xy / (2 * rBX * c2xy),
      -xy / (2 * rAY * c2xy),
      -xy / (2 * rBY * c2xy),
       xy / (rAB * c2xy)
    )

    idx6 <- list(c(ix, iy), c(ia, ix), c(ib, ix),
                 c(ia, iy), c(ib, iy), c(ia, ib))

    Sig <- matrix(0, 6, 6)
    for (m in seq_len(6))
      for (n in m:6) {
        v <- .cov_r(R, idx6[[m]][1], idx6[[m]][2],
                       idx6[[n]][1], idx6[[n]][2], N)
        Sig[m, n] <- Sig[n, m] <- v
      }

    var_diff <- as.numeric(t(g) %*% Sig %*% g)
    if (var_diff <= 0) { pvals[s] <- NA; next }

    pvals[s] <- 2 * pnorm(-abs(dz / sqrt(var_diff)))
  }

  pvals
}
