prob_Rr <- function(R, N, no) {
  # no, le nombre de dimensions exclues est 1 pour les corrélations résiduelles
  # où les résidus sont hors de l'hyper dagonale et hors du facteur exclu
  dl  <- N - 1 - no
  RR  <- 1 - R * R
  res <- rep(NA_real_, length(R))

  manq <- is.na(RR)
  neg  <- !manq & RR < 0
  ok   <- !manq & RR >= 0

  if (any(manq)) {
    warning("prob_Rr: ", sum(manq), " corrélation(s) NA reçue(s) — p-valeur(s) fixée(s) à NA.")
browser()
      }

  res[ok] <- 2 * pt(-abs(R[ok] * sqrt(dl / RR[ok])), dl)
  if (any(neg)) {
    cat("prob_Rr: |r| > 1 pour", sum(neg), "valeur(s) :",
        round(R[neg], 3), "\n")
    res[neg] <- 99.999  # artefact d'estimation — non significatif mais visible
  }

  # res reste NA là où R est NA
  return(res)
}
