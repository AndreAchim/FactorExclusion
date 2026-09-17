asGrappes <- function(AS) {
  # Calcule les distances (en X2)
  # Effectue les regroupements (linkage)
  # Établit le seuil de coupure en termes de probabilité (prob d'erreur de type I)
  # du lien plutôt que celle du pire
  # Détermine quels liens sont dans la fourchette AS$seuils[1]=.001 à AS$seuils[2]=.25
  # et crée des branches en fonction de cela, enlevant les liaisons dans l'ordre inverse de leur formation
  AS <- asDistances(AS)
  Z <- linkage(AS$Dist)
  nv <- length(AS$pertinent)
  Gr <- vector("list", nv)
  for (k in seq_len(nv)) {
    Gr[[k]] <- AS$pertinent[k]
  }
  nZ <- nrow(Z)
  pr <- numeric(nZ)
  nz <- 0
  for (k in seq_len(nZ)) {
    a <- Gr[[Z[k, 1]]]
    b <- Gr[[Z[k, 2]]]
    nx <- length(a) * length(b)
    nz <- nz + 1
    #    pr[k] <- 1 - pchisq(Z[k,3], nv - 2)^nx  # probabilité ajustée
    # la distance est maintenant -log10(pa) où pa est déjà la probabilité ajustée
    # Cette probabilité ne serait plus requise mais on la laisse pour compatibilité historique avec la suite
    #    p <- 10^-Z[k,3]
    p <- probCorr(Z[k,3],AS$N)$p
    pr[k] <- 1 - (1-p)^nx
    Gr[[nv + nz]] <- sort(c(a, b))
  }
  AS$GrBrut <- Gr[(nv + 1):length(Gr)]
  AS$Z <- cbind(Z, pr)
  if (max(pr) < AS$seuils[1]) {
    AS$VG[[1]] <- list(Gr = NULL, Creat = 'Initial (vide)', Parent = 0, coplan = NULL)
    return(AS)
  }
  if (any(pr > AS$seuils[2]))
  c <- tail(which(pr > AS$seuils[2]), 1)
  else
    c <- -nv
  Gr0 <- Gr[seq_len(nv + c)]
  res <- asNettoie(Gr0)
  
  #POC: Here is start to split into two according to seuils ####
  
  AS$VG[[1]] <- list(Gr = res$Gr, reste = res$reste, Creat = if (nv > 0) 'Initial' else 'Initial (vide)', Parent = 0)
  f <- which(pr > AS$seuils[1] & pr < AS$seuils[2])
  
  if (length(f) != 0) {
    
    ff <- subsets(f)
    ff <- exclutSequentiels(Z, ff)
    
    ng <- length(AS$VG[[1]]$Gr)
    mg <- 1
    for (k in seq_along(ff)) {
      AS$VG[[length(AS$VG) + 1]] <- AS$VG[[1]]
      AS$VG[[length(AS$VG)]]$LienInclus <- Gr[ff[[k]] + nv]
      AS$VG[[length(AS$VG)]]$Creat <- 'Grappes'
      AS$VG[[length(AS$VG)]]$Parent <- 1
      Gr1 <- c(Gr0, Gr[ff[[k]] + nv])
      if (length(AS$VG[[length(AS$VG)]]$Gr) < AS$minFct && c > 0) {
        AS$VG <- AS$VG[-length(AS$VG)]
      } else {
        res <- asNettoie(Gr1)
        AS$VG[[length(AS$VG)]]$Gr <- res$Gr
        AS$VG[[length(AS$VG)]]$reste <- res$reste
      }
      if (length(AS$VG[[length(AS$VG)]]$Gr) > ng) {
        ng <- length(AS$VG[[length(AS$VG)]]$Gr)
        mg <- length(AS$VG)
      }
    }
    if (mg > 1) {
      AS$VG[[mg]]$Parent <- 2
    }
    
    o <- c(mg, setdiff(seq_along(AS$VG), mg))
    AS$VG <- AS$VG[o]
  }
  return(AS)
}
