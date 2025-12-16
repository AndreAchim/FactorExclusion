logit_params <- function(rAtt,noms){  # vecteur des corrélations partielles moyennes attendues, Ar ou Br
  # nomsA ou nomsB
  r <- rep(abs(rAtt),200)
  parm <- list()
  k <- 0
  for (nom in noms){
    dt <- get(nom)[[2]][,13:18,]
    moyProb <- apply(dt,3,rowMeans)
    d <- as.vector(moyProb)<.05
    rd <- data.frame(d=d,r=r)
    fit <- glm(d~r,rd,family=binomial("logit"))
    fit=summary(fit)$coefficients
    fit <- cbind(fit,c(min(abs(rAtt)),max(rAtt)))
    colnames(fit) <- c("B","se","z","p","MiMa")
    parm <- append(parm,list(fit=fit))
  }
  names(parm) <- noms
  return(parm)
}