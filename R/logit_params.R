logit_params <- function(rAtt,noms){  # vecteur des corrélations partielles moyennes attendues, Ar ou Br
  # nomsA ou nomsB
  # refait le 2 janvier 2026
  # rAtt est maintenant une matrice
  Rm <- .6
  SF <- get("ST_Fl")
  paire <- SF$paire[[2]]
  np <- dim(paire)[2]
  R <- rep(NA,np)
  rAtt <- abs(rAtt)
  for (k in 1:np)
    R[k] <- rAtt[paire[1,k]-2,paire[2,k]-2]  # v1 et v2 n'existent pas dans rAtt
  r <- rep(R,1000)
  minR <- min(r)
  maxR <- max(r)
  parm <- list()
  k <- 0
  for (nom in noms){
    # dt <- log(get(nom)[[2]][,c(16,18),])
    # moyProb <- exp(apply(dt,3,rowMeans))
    N <- as.numeric(substr(nom,2,nchar(nom)))
    isem <- sqrt((N-3)*(2+2*Rm))
    dt <- atanh(get(nom)[[2]][,c(10,12),])
    moyProb <- 2*pnorm(-isem*abs(apply(dt,3,rowSums)))
    d <- as.vector(moyProb)<.05
    rd <- data.frame(d=d,r=r)
    fit <- glm(d~r,rd,family=binomial("logit"))
    fit=summary(fit)$coefficients
    fit <- cbind(fit,c(minR,maxR))
    colnames(fit) <- c("B","se","z","p","MiMa") # fir à 2 lignes, une pour chaque paramètre
    parm <- append(parm,list(fit=fit))
  }
  names(parm) <- noms
  return(parm)
}