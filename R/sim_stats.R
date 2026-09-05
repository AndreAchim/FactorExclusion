sim_stats <- function(fe, seuil = .05, N = NULL) {
  # Calcule les statistiques descriptives sur les simulations produites par FE_methodes2.
  #
  # Arguments :
  #   fe    : liste produite par FE_methodes2
  #   seuil : seuil de significativité (défaut .05)
  #   N     : taille d'échantillon (si NULL, utilise fe$N)
  #
  # Retourne : list(CR, sgCR, nOK)
  #   CR    : matrice 6 × np — statistiques de biais sur atanh des corrélations résiduelles
  #   sgCR  : matrice 2 × np — taux de détection à alpha = .05 et .01
  #   nOK   : nombre de simulations sans condition de Heywood
  if (is.null(N)) N <- fe$N
  OK  <- which(fe$Heywood == FALSE)
  nOK <- length(OK)

  fe_stats <- function(dt, att, pr) {
    att <- atanh(att)
    np  <- length(att)
    st  <- matrix(nrow = 6, ncol = np)
    sg  <- matrix(nrow = 2, ncol = np)
    for (k in seq_len(np)) {
      adt   <- atanh(dt[, k]) - att[k]
      tt    <- t.test(adt)
      st[, k] <- c(tt$statistic, tt$parameter, tt$p.value,
                   tanh(c(tt$estimate, tt$conf.int[1:2])))
      sg[1, k] <- mean(pr[, k] < .05, na.rm = TRUE)
      sg[2, k] <- mean(pr[, k] < .01, na.rm = TRUE)
    }
    return(list(st = st, sg = sg))
  }

  aCR <- fe_stats(fe$CR[OK, ], fe$CRatt, fe$pCR[OK, ])
  return(list(CR = aCR$st, sgCR = aCR$sg, nOK = nOK))
}
