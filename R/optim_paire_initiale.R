optim_paire_initiale <- function(AS, ab) {
  # Optimise le poids p du contraste -a + p*b pour minimiser R².
  autres <- setdiff(AS$pertinent, ab)
  col    <- AS$GS[, ab]
  TEMOINS <- AS$GS[, autres]
  Rinv   <- solve(t(TEMOINS) %*% TEMOINS)
  nt     <- length(autres)
  aa     <- rep(NA, nt)
  bb     <- rep(NA, nt)
  rAB    <- AS$R[ab[1], ab[2]]
  for (k in 1:nt) {
    if (rAB * prod(AS$R[ab, k]) <= 0) {
      aa[k] <- 1
      bb[k] <- 0
    } else {
      aa[k] <- rAB * AS$R[ab[1], k] / AS$R[ab[2], k]
      bb[k] <- rAB * AS$R[ab[2], k] / AS$R[ab[1], k]
      if (is.nan(aa[k]) || is.nan(bb[k]))
        stop("NaN dans le calcul des poids initiaux.")
    }
  }
  aa <- sign(rAB) * sqrt(aa)
  bb <- sqrt(bb)
  pp <- median(bb / aa)
  if (is.nan(pp)) stop("Le poids initial médian est NaN.")
  out      <- optim(pp, crit_R2, gr = NULL, col, TEMOINS, Rinv, method = "BFGS")
  R2       <- out$value
  prob     <- prob_R2(R2, nt, AS$N)
  po       <- out$par
  contrast <- col %*% c(-1, po)
  corr     <- t(sc1(contrast)) %*% TEMOINS
  return(list(crit = R2, po = po, prob = prob, corr = corr, contrast = contrast))
}
