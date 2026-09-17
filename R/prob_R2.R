prob_R2 <- function(R2,k,N){
  dl <- N-k-1
  Fr <- dl*R2/(k*(1-R2))
  return(1-pf(Fr,k,dl))
}