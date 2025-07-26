rapport_z2 <- function(Z2,aa=NULL){
# z2 est (npaires,nrep), aa est la liste des rangs de paires avec annulation attendue
  zcrit <- qchisq(c(.75,.999),1)
  paires <- Z2$paires
  z2 <- Z2$z2
#  browser()
  for (p in 1:nrow(z2)){
    ni <- which(is.infinite(z2[p,]))
    if (length(ni)>0)
      z2[p,ni] <- NaN
    nb[5] <- 0
    nb[4] <- length(z2[p,])
    nb[3] <- sum(!is.nan(z2[p,]))
    nb[2] <- sum(z2[p,] < zcrit[2],na.rm=TRUE)
    nb[1] <- sum(z2[p,] < zcrit[1],na.rm=TRUE)
    nb[2:5] <- diff(nb)
    if (length(ni)>0)
      z2[p,ni] <- Inf
    nb[4] <- sum(is.nan(z2[p,]))
    nb[5] <- sum(is.infinite(z2[p,]))
    if (any(aa==p)){
      mima <- max
      mim <- "max"
      z2[p,z2[p,]==Inf] <- -Inf
    } else {
      mima <- min
      mim <- "min"
    }
    legen <-" Nb(>.25 .25:.001 <.001 NaN r(ns)): "
    sp <- sprintf("%3.0f,%3.0f, %s = %3.4f %s  %3d %3d %3d %3d %3d",
                  paires[p,1],paires[p,2],
                  mim,mima(z2[p,],na.rm=TRUE),legen,nb[1],nb[2],nb[3],nb[4],nb[5])
    print(sp,quote=FALSE)
  }
}
