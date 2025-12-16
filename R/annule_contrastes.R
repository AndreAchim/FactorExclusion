annule_contrastes <- function(AS){
  if (length(AS$excl_co)<3)
    stop("Le champ %excl_co n'a pas trois éléments")
  nc <- ncol(AS$excl_co[[1]])
  np <- nc * (nc-1) / 2
  paire <- matrix(0,np,2)
  crit <- rep(0,np)
  poids <- rep(0,np)
  proj <- matrix(0,np,nc-2)
  p <- 0
  for (j in 1:(nc-1))
    for (k in (j+1):nc){
      p <- p+1
      paire[p,] <- c(j,k)
      out <- exclude_fact(AS$excl_co[[1]][,j,drop=FALSE],AS$excl_co[[2]][,k,drop=FALSE],AS$excl_co[[3]][,-c(j,k)])
      crit[p] <- out$crit * AS$N
      poids[p] <- out$weights
      proj[p,] <- out$proj
    }
  return(list(paire=paire,crit=crit,poids=poids,proj=proj))
}  