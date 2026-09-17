rapport_exclusion <- function(AS){
  # ncar, si présent est le nombre de caractères dans les titres des colonnes de AS$GS
  # Toutes et seules les colonnes ainsi sélectionnées sont traitées par paire de variable
  # traitée (identifiée par la première lettre du nom) mais seulement les combinaisons qui
  # ne partagent pas de bruit (de lettres dans leur noms)
  ra <- r_atanh(AS)
  dra <- ra-rowMeans(ra)
  dra <- SC1(t(dra))
  RR <- t(dra) %*% dra
  tRR <- sum(RR)
  iSeTot <- sqrt((AS$N-3)/tRR)  # inverse d'erreur type du total
  kRR <- nrow(RR)
  Ratanh <- (tRR-kRR)/(kRR*(kRR-1)) # estimation de corrélations entre les atanh des corrélations partielles
  # browser()
  noms <- colnames(AS$GS)
  ncar <- nchar(noms[length(noms)])
  cl <- which(nchar(noms)==ncar)   # rangs des colonnes à ncar caractères
  noms <- noms[cl]
  deb <- substring(noms,1,1)
  var <- unique(deb)
  nv <- length(var)
  CR <- matrix(NA,nrow=nv,ncol=nv)
  R <- matrix(NA,nrow=nv,ncol=nv)
  P <- matrix(NA,nrow=nv,ncol=nv)
  RR <- array(NA,c(nv,nv,4))
  NOMS16 <- list()
  CR16 <- NULL
  R16 <- NULL
  NN <- 0
  for (v1 in 1:(nv-1)){
    c1 <- AS$GS[,cl[deb==var[v1]]]
    xnoms <- noms[deb==var[v1]]
    noms1 <- strsplit(xnoms,"")
    for (v2 in (v1+1):nv){
      c2 <- AS$GS[,cl[deb==var[v2]]]
      ynoms <- noms[deb==var[v2]]
      noms2 <- strsplit(ynoms,"")
      cr <- NULL
      r <- NULL
      paire <- NULL
      for (k1 in 1:length(noms1)) {
        for (k2 in 1:length(noms2)){
          qq <- length(intersect(noms1[[k1]],noms2[[k2]]))
          if (qq==0){
            paire <- paste(paire,paste0(xnoms[[k1]],ynoms[[k2]]))
            cr <- c(cr,t(c1[,k1]) %*% c2[,k2])   # covariance résiduelle par paire sans bruit partagé
            r <- c(r,sc1(t(c1[,k1])) %*% sc1(c2[,k2]))  # corrélation partielle correspondante
          }
        }
      }
      if (!is.null(paire)){
        NN <- NN+1
        NOMS16[[NN]] <- paire
        CR16 <- cbind(CR16,cr)
        R16 <- cbind(R16,r)
        CR[v2,v1] <- mean(cr)
        R[v1,v2] <- mean(r)
        R[v2,v1] <- tanh(mean(atanh(r)))
        # CR[v1,v2] <- 2*pnorm(-abs(kRR*R[v2,v1]*iSeTot))
        CR[v1,v2] <- 2*pnorm(-abs(sum(atanh(cr))*iSeTot))
      }
    }
  }
  out <- nommer(list(CR=CR,R=R),var)
  out$Rm <- Ratanh
  out$NOMS16 <- NOMS16
  out$CR16 <- CR16
  out$R16 <- R16
  out$Rm <- Ratanh
  return(out)
}

nommer <- function(MAT,noms){
  for (m in 1:length(MAT)){
    rownames(MAT[[m]]) <- noms
    colnames(MAT[[m]]) <- noms
  }
  return(MAT)
}