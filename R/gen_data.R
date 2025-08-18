library(MASS)
library(expm)
gen_data <- function(F,N,R=NULL,with_var=FALSE){
# F est matrice de patrons de nv rangées
# R, si présente, est la matrice de corrélation des facteurs
# si with_var est vraie, la dernière colonne de F est la variance des unicités
# autrement la variance d'unicité est 1 moins la somme des carrés sur les colonnes
  if (!is.matrix(F)){
    nf <- 1
    nv <- length(F)
    unic=1-F*F
    f <- as.matrix(F)
  } else {
    nv <- nrow(F)
    nf <- ncol(F)
    if (with_var){
      unic <- F[,nf]
      F <- F[,-nf]
      nf <- nf-1
    } else
      unic <- 1-rowSums(F * F)
  }
  if (is.null(R))
    R <- diag(rep(1,nf))
  src_cov <- diag(nf+nv)
  src_cov[1:nf,1:nf] <- R
  G <- cbind(F %*% sqrtm(R),diag(sqrt(unic)))
  dat <- MASS::mvrnorm(N,rep(0,nf+nv),src_cov) %*% t(G)
}
  