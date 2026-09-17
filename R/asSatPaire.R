# Cpaires, Ppaires, R
asSatPaire <- function(AS, v){
  inverse <- (v[1] > v[2]) # Techniquement impossible ???
  if(inverse) v <- v[c(2, 1)]
  
  r <- AS$R[v[1], v[2]]
  pa <- which(apply(AS$Cpaires == matrix(v, ncol = ncol(AS$Cpaires), 2, 
                                         byrow = FALSE), 2, all))
  p <- AS$Ppaires[pa]
  #crit <- AS$Crit[pa]
  if((p * r) < 0){
    #sat <- c(0, 0)
    DD <- warnings();if(length(DD)>0) browser()
    sat <- sqrt(abs(r)) * c(1, sign(r))
  } else {
    s1 <- sqrt(r/p)
    s2 <- r/s1
    if(inverse) {
      sat <- c(s2, s1)
    } else{
      sat <- c(s1, s2)
    }
  }
  return(sat)
  # POC : should return crit and AS??
}
