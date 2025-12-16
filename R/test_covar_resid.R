test_covar_resid <- function(AS,cx,cy,v_excl=c(1,2),srce=NULL){
  # cx et cy sont des rangs de variables; v_excl donne la liste des covariables d'exclusion de dimension. 
  # Forme pour chaque cas les contrastes AS$excl_weig[v_excl[1]][c(cx,cy)] et AS$excl_weig[v_excl[2]][c(cx,cy)]
  # Moyenne pour chaque cas les produits (cx,cy) de listes différentes
  # srce n'est pas utilisé mais pourra être implanté quand on voudra annuler un facteur sur des contrastes
  # bien que ceci pourrait être fait, j'imagine, en mettant ces contrastes dans AS$zdat
  # Si L est la longueur de v_excl, produit L*(L-1) produits croisés pour chaque sujet, en confectionne la
  # matrice de leurs produits croisés moyens (sans exclusion des moyennes), en extrait la première composante
  # principale, porte la somme de ses poids à 1 et l'applique aux L(L-1) produits de chaque sujet
  # comme scores de chaque sujet sur cette composante. Applique le t de Student pour tester
  # si la moyenne diffère de 0. Cette moyenne serait la covariance résiduelle.
  if (is.null(dim(AS$zdat)))
    AS <- cree_zdat(AS)
  x <- which(cx==AS$excl_cibles)
  y <- which(cy==AS$excl_cibles)
  if (is.null(x) || is.null(y)) error("cx et cy doivent être des rangs de variables inclus dans AS$util")
  ncov <- length(v_excl)
  nc1 <- ncov+1
  X=matrix(NA,nrow=AS$N,ncol=nc1)
  X[,nc1] <- AS$zdat[,cx]
  Y=matrix(NA,nrow=AS$N,ncol=nc1)
  Y[,nc1] <- AS$zdat[,cy]
  for (k in 1:ncov) {
    X[,k] <-  AS$zdat[,cx] - AS$excl_weig[[k]][x]*AS$zdat[,AS$excl_var[k]]
    Y[,k] <-  AS$zdat[,cy] - AS$excl_weig[[k]][y]*AS$zdat[,AS$excl_var[k]]
  }
  PROD <- matrix(NA,nrow=AS$N,ncol=ncov*(ncov+1))
  p <- 0
  for (i in 1:nc1){
    for (j in 1:nc1){
      if (i != j){
        p <- p+1
        PROD[,p] <- X[,i] * Y[,j]
      }
    }
  }
  moy <- colMeans(PROD)
  if (any(sign(moy)!=sign(moy[1])))
    return(.5)
#  CR <- svd(PROD,nu=1,nv=1)
#  tt <- t.test(CR$u)
  tt <- t.test(apply(PROD[,3:4],1,mean))
  return(tt$p.value)
#  return(list(CR=CR,ttest=tt,moyPROD=moy))
#  browser()
#   po <- do.call(rbind,AS$excl_weig)[,c(x,y)]
# #  dl <- AS$N-2
# 
# #  stat <- list()
#   s <- 0
#   n <- ns*(ns-1)/2
#   val <- matrix(NA,nrow=AS$nv,ncol=6)
#   val[,6] <- AS$GS[,cy]
#   val[,5] <- AS$GS[,cx]
#   #  paire <- matrix(c(6,1,5,2,4,1,3,2,6,3,5,4),nrow=2)
#   p1 <- c(6,5,4,3,6,5)
#   p2 <- c(1,2,1,2,3,4)
#   CR <- matrix(NA,nrow=6,ncol=n)
#   #  pond <- matrix(NA,nrow=2,ncol=n)
#   r <- matrix(NA,nrow=6,ncol=n)
#   pr <- matrix(NA,nrow=6,ncol=n)
#   #browser()  
#   for (j in 1:(ns-1)){
#     val[,1] <-AS$excl_cont[[AS$excl_var[j]]][,x]
#     val[,2] <-AS$excl_cont[[AS$excl_var[j]]][,y]
#     for (k in (j+1):ns){
#       val[,3] <-AS$excl_cont[[AS$excl_var[k]]][,x]
#       val[,4] <-AS$excl_cont[[AS$excl_var[k]]][,y]
#       # pour chaque paire de variable d'exclusion
#       s <- s+1
#       # v1 <- j
#       # v2 <- k
#       for (v in 1:6) {
#         # C1x <-  AS$zdat[,cx] - AS$excl_weig[[1]][x]*AS$zdat[,AS$excl_var[1]]
#         # C1y <-  AS$zdat[,cy] - AS$excl_weig[[1]][y]*AS$zdat[,AS$excl_var[1]]
#         # C2x <-  AS$zdat[,cx] - AS$excl_weig[[2]][x]*AS$zdat[,AS$excl_var[2]]
#         # C2y <-  AS$zdat[,cy] - AS$excl_weig[[2]][y]*AS$zdat[,AS$excl_var[2]]
#         # pour chaque appariment d'une cible avec une variable d'exclusion
#         # Cx <- AS$excl_cont[[AS$excl_var[v1]]][,x]
#         # Cy <- AS$excl_cont[[AS$excl_var[v2]]][,y]
#         CR[v,s] <- t(val[,p1[v]]) %*% val[,p2[v]]
#         r[v,s] <- t(sc1(val[,p1[v]])) %*% sc1(val[,p2[v]])
#         tr <- r[v,s]*sqrt(dl/(1-r[v,s]*r[v,s]))
#         pr[v,s] <- 2*pt(-abs(tr),dl)
#         # v1 <- k
#         # v2 <- j
#       }
#       # # calculer la pondération pour le produit de cette paire de contrastes
#       # p1 <- AS$excl_weig[[v1]]*sat[v1]
#       # p2 <- AS$excl_weig[[v2]]*sat[v2]
#       # kk=abs(AS$R[cx,cy] - p1[x]*p2[y])
#       # pond[1,s] <- 1/sqrt(1-p1[x]^2-kk+(AS$excl_weig[[v1]][x]^2*(1-sat[v1]^2)))
#       # pond[2,s] <- 1/sqrt(1-p2[y]^2-kk+(AS$excl_weig[[v2]][y]^2*(1-sat[v2]^2)))
#       # CR_ <- sum(CR*pond)/sum(pond)
#       # if (all(sign(CR)==sign(CR[1,1]))) {
#       #   pr_ <- sum(pr*pond)/sum(pond)
#       #   PR <- 2*pnorm(-abs(sum(qnorm(pr))/sqrt(ns)))
#       # } else
#       if (any(sign(CR)!=sign(CR[1,1]))){
#         pr <- rep(.5,max(6,n))
#       }
#     }
#   }
#   #  return(list(CR_=CR_,CR=CR,r=r,PR=PR,pr_=pr_,pr=pr,pond=pond))
# #   return(list(CR=CR,r=r,pr=pr))
  # return(c(CR,r,pr,AS$satur_var))
  #  return(c(mean(CR),mean(r),mean(pr)))
}
