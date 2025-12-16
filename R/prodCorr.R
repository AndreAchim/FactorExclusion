prodCorr <- function(r,N,corr_dl=0){
  dl <- N-2-corr_dl
  t <- r*(sqrt(dl/(1-r*r)))
  p <- 2*pt(-abs(t),dl)
  return(list(p=p,t=t,dl=dl))
}