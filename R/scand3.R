scand3 <- function(R,N=NULL){
  # if (!isSymmetric(R)){
  #   N=nrow(R)
  #   R <- cor(R)
  # }
  AS <- init_SCA(R,N)
  out <- list()
  encore <- TRUE
  k <- 0
  vDir <- 0
  while (encore) {
    k <- k+1
    tst <- seq_k_dim(AS,vDir)
    if (k==1) vDir <- tst$meil
    out[[length(out)+1]] <- tst
    prb <- tst$out[,ncol(tst$out)]
    # print(sum(prb<.01))
    encore <- sum(prb<.01)>=3
    # encore <- tst$prob<.05
    vDir <- c(tst$prochain,vDir)
  }
  return(out)
}