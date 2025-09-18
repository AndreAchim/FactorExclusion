cree_zdat <- function(AS){
  if (is.null(AS$zdat)){
    AS$zdat <- AS$dat
    di <- dim(AS$dat)
    if (di[1] > di[2]){ 
      for (k in 1:di[2]){
        zd <- AS$zdat[,k]
        zd <- zd - mean(zd)
        s <- 1/ sqrt(t(zd) %*% zd)
        AS$zdat[,k] <- zd * as.vector(s)
      }
    }
  }
  return(AS)
}