dist_z2 <- function(F,N=75,nrep=200){
# inspection des distributions de z2 et en particulier du critère d'optimisation
# le z2 maximum
# nrep <- 200
# N <- 75
# F <- matrix(c(.6,.4,.8,.4,.7,0,0, 0,0,0,.7,.4,.5,.8),ncol=2)
  nv <- nrow(F)
  nf <- ncol(F)
  paires <- paires_bifactorielles(F)
  paires <- rbind(c(1,2,1),paires) # ajouter une paire à annulation correcte possible
  np <- nrow(paires)
  di <- c(nv-2,np,nrep)
  z2 <- array(rep(NaN,prod(di)),di)
  for (k in 1:nrep){
    for (j in 1:np){
      a <- paires[j,1]
      b <- paires[j,2]
#      cible=setdiff(1:nv,c(a,b))
      rdat <- chol(cor(gen_data(F,N)))
#      limi <- init_lim(rdat,a,b)
#      browser()
#      p0 <- mean(limi$limi[1])
#      print(c(k,j,projections(p0,rdat,a,b),p0))
#      if (max(limi$limi)>4)
#        limites <- c(-200,200)
#        else
#          limites <- c(p0-.3,p0+.3)
      out <- optimize(max_proj,limites,rdat,a,b,fct=max_abs)
      #       suppressWarnings({
#       out <- optim(1, fn = pair_cancel, G = rdat, method = "Nelder-Mead",
#                       combine = c(a,b),
#                       cible = cible,
#                       control = list(maxit = 1e9,
#                                      factr = 1e-9,
#                                      pgtol = 1e-9))
#       })
# #      Po[, e] <- fopt$par
#      pr <- asCrit(out$par,rdat,c(a,b), cible, to.opt = FALSE)$cor[-c(a,b)]
      pr <- projections(out$minimum,rdat,a,b)
      print(c(k,j,pr))
      if (abs(max(pr)+min(pr))>.1) browser()
      z2[,j,k] <- (N-1) * pr
    }
  }
  return(list(z2=z2,paires=paires))
}