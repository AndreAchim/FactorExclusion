FE_methodes2 <- function(F,N,ote=1:2,nrep=1000){
  a <- ote[1]
  b <- ote[2]
  nv <- nrow(F)
  nf <- ncol(F)
  if (N<0){
    N <- -N
    set.seed(N)
  }
  if (is.null(rownames(F)))
    rownames(F) <- letters[1:nv]
  
  # calculer les corrélations résiduelles attendues
  G <- as.matrix(F[-ote,2:nf])
  R <- G %*% t(G)
  lR <- lower.tri(R)
  uR <- upper.tri(R)
  CRatt <- R[lR]
  
  # réserver l'espace
  ng <- nv - length(ote)
  nc <- ng*(ng-1)/2
  CR <- matrix(nrow = nrep,ncol=nc)  # corrélations résiduelles
  Heywood <- logical(nrep)
  pCR <- matrix(nrow = nrep,ncol=nc)
  lR0 <- which(abs(CRatt)<1e-5)
  # lR0 <- which(CRatt[a0])
  n0 <- length(lR0)
  # CP <- matrix(nrow = nrep,ncol=n0)  # corrélations partielles
  # pCP <- matrix(nrow = nrep,ncol=n0)
  # CR0 <- matrix(nrow = nrep,ncol=nc)  # corrélations résiduelles par annulation
  # pCR0 <- matrix(nrow = nrep,ncol=nc)
  # dx <- matrix(nrow = nrep, ncol=2)
  # nFP <- matrix(0,nrow=3,ncol=5)
  # zro <- c(ote,ote+nv)
  remplaces <- 0
  for (k in 1:nrep){
    refait <- TRUE
    while(refait){
      dat <- gen_data(F,N)$dt
      R <- cor(dat)
      AS <- pre_partition(R,N,ote,seuil=.05)
      refait <- is.null(AS)
      if (refait)  remplaces <- remplaces + 1
    }
    # dx[k,1] <- R[a,b]
    # if (n0){  # seulement s'il y a des corrélations résiduelles attendues nulles
    #   Rp <- corr_partielles(R,N,ote)$Rp
    #   lRp <- Rp[lower.tri(Rp)]
    #   CP[k,] <- lRp[lR0]
    #   tRp <- t(Rp)
    #   lRp <- tRp[lower.tri(tRp)]
    #   pCP[k,] <- lRp[lR0]
    # }
    pR <- partitionne_R(AS,ote)
    Rp <- pR$corr_resid
    # dx[k,2] <- ifelse(pR$sat$Heywood[seuil],1,0)
    lRp <- Rp[lower.tri(Rp)]
    CR[k,] <- lRp  #[lR0]
    tRp <- t(Rp)
    pCR[k,] <-  tRp[lower.tri(tRp)]
    Heywood[k] <- pR$Heywood
    if ((k %% 100)==0) cat(k,"/",nrep," rempl:",remplaces,"\n")
  }
  # sig <- t(matrix(rowSums(0+(out<.05)),ncol=2))
  # Nb <- rowSums(0+(out[1:nc,]<5))
  # IC99 <- floor(.05*Nb+2.56*sqrt(.95*.05*Nb))
  # sig <- rbind(sig,Nb,IC99)
  # return(list(sig=sig,out=out))
  # return(list(CR=CR,pCR=pCR,CP=CP,pCP=pCP,CR0=CR0,pCR0=pCR0,CRatt=CRatt,dx=dx))
  # return(list(CR=CR,pCR=pCR,CP=CP,pCP=pCP,CRatt=CRatt,N=N,dx=dx))
  return(list(CR=CR,pCR=pCR,CRatt=CRatt,N=N,Heywood=Heywood))
}