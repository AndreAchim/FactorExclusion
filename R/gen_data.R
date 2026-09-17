library(MASS)
library(expm)
gen_data <- function(F,N,R=NULL,with_var=FALSE,ON=FALSE){
  # F est matrice de patrons de nv rangées
  # R, si présente, est la matrice de corrélation des facteurs
  # si with_var est vraie, la dernière colonne de F est la variance des unicités
  # autrement la variance d'unicité est 1 moins la somme des carrés sur les colonnes
  if (!is.matrix(F)){  # si une seule colonne à F
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
    } else {
      unic <- 1-rowSums(F * F)
    }
  }
  
  if (!is.null(R)){
    if (isSymmetric(R))
      R <- t(chol(R))
    else if (sum(abs(R(lower.tri(R)))==0))
      r <- t(R)
    F <- F %*% R
  }
  src_cov <- diag(nf+nv)
  G <- cbind(F,diag(sqrt(unic)))
  srce <- MASS::mvrnorm(N,rep(0,nf+nv),src_cov,empirical=ON)
  # srce <- matrix(rnorm(N*(nf+nv)),nrow=N)
  dat <- srce %*% t(G)
  colnames(dat) <- letters[1:nv]
  return(list(dt=dat,srce=srce))
}
