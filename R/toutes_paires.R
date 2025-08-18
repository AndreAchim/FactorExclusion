toutes_paires <- function(F){
# F est de préférence une matrice de patrons plutôt que de structure
# Cela permet de retourner aussi la liste des paires de variables exclusives à un même facteur
# Celles-ci sont signalées par le rang de la première variable de la paire rapportéen négatif
  nv <- nrow(F)
  nf <- ncol(F)
  np <- nv*(nv-1)/2
  paires <- matrix(0,nrow=np,ncol=2)
  p <- 0
  for (v1 in 1:(nv-1))
    for (v2 in (v1+1):nv){
      p <- p+1
      paires[p,] <- c(v1,v2)
      w1 <- which(F[v1,] != 0)
      if (length(w1) == 1){
        w2 <- which(F[v2,] != 0)
        if (length(w2) ==1 && w1==w2)
          paires[p,1] <- -v1
      }
    }
  return(paires)
  # browser()
  # return(paires[1:p,])
}
