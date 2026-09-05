as_annule_paire <- function(AS, k) {
  # Optimise les poids de la paire k (pointée dans AS$Cpaires) pour minimiser
  # le signal résiduel. Remplit AS$Crit[k], AS$Ppaires[k], AS$Prob[k],
  # AS$satPaires[,k].
  melange <- AS$Cpaires[, k]
  pd      <- probNonDoublet(AS, melange)
  r       <- AS$R[melange[1], melange[2]]
  if (pd < 0.05) {
    ap <- optim_paire_initiale(AS, melange)
  } else {
    # La paire semble dépendre d'un facteur doublet
    AS$doublet <- rbind(AS$doublet, melange)
    out  <- probCorr(r, AS$N)
    cont <- AS$GS[, melange] %*% c(-1, 1)
    corr <- t(sc1(cont)) %*% sc1(AS$GS[, -melange])
    cible <- setdiff(AS$pertinent, melange)
    R2   <- crit_R2(1, AS$GS[, melange], AS$GS[, cible])
    prob <- prob_R2(R2, length(cible), AS$N)
    ap   <- list(prob = prob, po = sign(r), crit = R2, corr = corr)
    ap$po <- sign(r)
  }
  if (abs(log(abs(ap$po))) > 5) {
    AS$Prob[k] <- .1^abs(log(abs(ap$po)))
  } else {
    AS$Prob[k] <- ap$prob
  }
  AS$Crit[k]    <- ap$crit
  AS$Ppaires[k] <- ap$po
  # Estimer les saturations
  if (r * ap$po > 0) {
    s1 <- sqrt(r * ap$po)
    s2 <- r / s1
  } else {
    s1 <- s2 <- NA
  }
  AS$satPaires[, k] <- c(s1, s2)
  return(AS)
}
