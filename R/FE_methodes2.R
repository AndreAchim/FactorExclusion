FE_methodes2 <- function(F, N, ote = 1:2, nrep = 1000) {
  # Simule nrep échantillons à partir de la matrice de saturations F
  # et calcule les corrélations résiduelles via la méthode FE.
  #
  # Arguments :
  #   F    : matrice de saturations (nv × nf)
  #   N    : taille d'échantillon (si négatif, fixe la graine = |N|)
  #   ote  : indices des deux variables ancres
  #   nrep : nombre de réplications
  #
  # Retourne : list(CR, pCR, CRatt, N, Heywood)
  a  <- ote[1]
  b  <- ote[2]
  nv <- nrow(F)
  nf <- ncol(F)

  if (N < 0) {
    N <- -N
    set.seed(N)
  }

  if (is.null(rownames(F)))
    rownames(F) <- letters[1:nv]

  # Corrélations résiduelles attendues (population)
  G      <- as.matrix(F[-ote, 2:nf])
  R_pop  <- G %*% t(G)
  lR     <- lower.tri(R_pop)
  CRatt  <- R_pop[lR]

  # Réserver l'espace
  ng  <- nv - length(ote)
  nc  <- ng * (ng - 1) / 2
  CR      <- matrix(nrow = nrep, ncol = nc)
  pCR     <- matrix(nrow = nrep, ncol = nc)
  Heywood <- logical(nrep)
  remplaces <- 0

  for (k in 1:nrep) {
    refait <- TRUE
    while (refait) {
      dat    <- gen_data(F, N)$dt
      R_obs  <- cor(dat)
      AS     <- pre_partition(R_obs, N, ote, seuil = .05)
      refait <- is.null(AS)
      if (refait) remplaces <- remplaces + 1
    }
    pR       <- partitionne_R(AS, ote)
    Rp       <- pR$corr_resid
    lRp      <- Rp[lower.tri(Rp)]
    CR[k, ]  <- lRp
    tRp      <- t(Rp)
    pCR[k, ] <- tRp[lower.tri(tRp)]
    Heywood[k] <- pR$Heywood
    if ((k %% 100) == 0)
      cat(k, "/", nrep, " rempl:", remplaces, "\n")
  }

  return(list(CR = CR, pCR = pCR, CRatt = CRatt, N = N, Heywood = Heywood))
}
