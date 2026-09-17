crit_R2 <- function(p,VAR,TEMOINS,Rinv=NULL){
  y <- sc1(VAR %*% c(-1,p))
  if (is.null(Rinv))
    Rinv <- solve(t(TEMOINS) %*% TEMOINS)
  z <- t(y) %*% TEMOINS
  R2 <- z %*% Rinv %*% t(z)
  return(R2)
}