corr_partielles <- function(R,N,ote){
  nv <- nrow(R)
  o <- c(setdiff(1:nv,ote),ote)
  no <- length(ote)
  S <- R[o,o]
  n1 <- nv-no
  un <- 1:n1
  de <- (n1+1):nv
  Sp <- S[un,un] - S[un,de] %*% solve(S[de,de]) %*% S[de,un]
  d <- diag(1/sqrt(diag(Sp)))
  Rp <- d %*% Sp %*% d
  diag(Rp) <- NA
  Rp[upper.tri(Rp)] <- prob_Rr(Rp[upper.tri(Rp)],N,no)
  return(list(Rp=Rp,Sp=Sp))
}