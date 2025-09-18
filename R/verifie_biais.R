verifie_biais <- function(N=200,nrep=100){
  F <- matrix(c(.8,.4,.8,.3,.6,0,0,.3,.8,.6),ncol=2)
  attendu <- t(matrix(c(F[3:5,1]/F[1,1],F[3:5,1]/F[2,1]),nrow=3))
  out <- array(dim=c(2,3,nrep))
  for (k in 1:nrep){
    dt <- gen_data(F,N)
    AS <- init_SCA(dt)
    AS <- cree_zdat(AS)
    AS <- factor_exclusion(AS,c(1,2))
    out[1,,k] <- AS$excl_weig[[1]]-attendu[1,]
    out[2,,k] <- AS$excl_weig[[2]]-attendu[2,]
  }
  return(list(att=attendu,poids=out,AS=AS))
}