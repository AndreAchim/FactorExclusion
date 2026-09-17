FE_sim <- function(F,N,ote=c(1,2),nrep=100){
  # F est une de patron(nv,nf)
  # N est la taille des échantillons et nrep le nombre de jeux de données simulées
  # ote désigne les deux variables d'exclusion
  # En sortie, on aura une liste de deux matrices (nrep,(nv-2)*(nv-3)/2)
  # $CR et $pCR pour les correlations résiduelles
  # et deux matrices $CP et $pCP (nrep,n0) où n0 est le nombre de
  # corrélations résiduelles attendues essantiellement nulles
  # if(length(warnings())>0) clear_warn()
  nv <- nrow(F)
  nf <- ncol(F)
  if (N<0){
    N <- -N
    set.seed(N*nf)
  }
  
  # calculer les corrélations résiduelles attendues
  G <- as.matrix(F[-ote,2:nf])
  R <- G %*% t(G)
  lR <- lower.tri(R)
  uR <- upper.tri(R)
  CRatt <- R[lR]
  
  # réserver l'espace
  nc <- (nv-2)*(nv-3)/2
  CR <- matrix(nrow = nrep,ncol=nc)
  pCR <- matrix(nrow = nrep,ncol=nc)
  a0 <- which(abs(CRatt)<1e-5)
  n0 <- length(a0)
  CP <- matrix(nrow = nrep,ncol=n0)
  
  # effectuer les nrep simulations
  # options(warn=1)
  for (rp in 1:nrep){
    out <- gen_data(F,N=N)
    dt <- out$dt
    AS <- init_SCA(dt,seuils=c(.6,.70)) # pour éviter les orphelines
    AS <- factor_exclusion(AS,ote)
    CC <- warnings();if(length(CC)>0) browser()
    RE <- rapport_exclusion(AS)
    Rm[rp] <- RE$Rm
    for (k in 1:2){   # attention si plus que 2
      result[[k]][,,rp] <- t(rbind(RE$CR16[,pa[[k]]],RE$R16[,pa[[k]]]))
    }
    #    result[[k]]
    # for (Li in SF$paires){
    #   #      a <- !a
    #   k <- k + 1
    #   for (p in 1:np[k]){
    #     out <- test_proj(AS,SF$paires[[k]][1,p],SF$paires[[k]][2,p],srce=srce)
    #     out <- c(out$CR,out$r,out$pr,AS$satur_var)
    #     #       if (a){
    #     #        aa <- test_cov_resid(AS,SF$paires[[k]][1,p],SF$paires[[k]][2,p])
    #     #          print(c(rp,k,p,aa,mean(out[13:18])))
    #     #        }
    #     # return(list(CR_=CR_,CR=CR,r=r,PR=PR,pr_=pr_,pr=pr,pond=pond))
    #     # result[[k]][p,,rp] <- c(out$CR_,mean(out$CR),mean(out$r),out$PR,out$pr_,mean(out$pr))
    #     #        if (length(out) != 18) browser()
    #     # browser()
    #     result[[k]][p,,rp] <- out
    #     # if (p==np[k] && all(sign(out[2:6])==sign(out[1]))){
    #     #   aa <-  exp(mean(log(out[13:18])))
    #     #   if (Li[1,1]==3 && aa<.05){
    #     #     V1 <- AS$zdat[,1]
    #     #     V2 <- AS$zdat[,2]
    #     #     V7 <- AS$zdat[,7]
    #     #     V8 <- AS$zdat[,8]
    #     #     X1 <- sc1(AS$zdat[,7]-AS$excl_w[[1]][5]*AS$zdat[,1])
    #     #     Y1 <- sc1(AS$zdat[,8]-AS$excl_w[[1]][6]*AS$zdat[,1])
    #     #     X2 <- sc1(AS$zdat[,7]-AS$excl_w[[2]][5]*AS$zdat[,2])
    #     #     Y2 <- sc1(AS$zdat[,8]-AS$excl_w[[2]][6]*AS$zdat[,2])
    #     #     AA <- t(srce[,1:3]) %*% cbind(V1,V2,V7,V8,X1,X2,Y1,Y2)
    #     #     #            print(sprintf("%3d %6.4f   %6.4f %6.4f %6.4f %6.4f",rp,aa,AA[1],AA[2],AA[3],AA[4]),quote=FALSE)
    #     #     print(sprintf("rp%d p%d Li%d_%d p%6.4f",rp,p,Li[1,p],Li[2,p],aa),quote=FALSE)
    #     #     print(cbind(t(AS$excl_weig[[1]][5:6]),t(AS$excl_weig[[2]][5:6]),t(AS$excl_crit[[1]][5:6]),t(AS$excl_crit[[2]][5:6])))
    #     #     #            if (aa<.01){
    #     #     print(AA)
    #     #     browser()
    #     #   }
    #     #    #           }
    #     # } else
    #     #   sdif[k] <- sdif[k]+1
    #     #        result[[k]] <- out
    #     #        if(rp==1) print(c(p,out))
    #     # }
    #     # browser()
    #     # result[[k,rp]] <- out
    #   }
    # }
    # #    AAA <- warnings();if(length(AAA)>0) browser()
  }
  return(list(result=result,Rm=Rm))
#  return(list(result=result,RES=RES))
}

prepare_poids_fixes <- function(AS,SF){
  v <- sort(unique(as.vector(SF$paires$zero))) +1  # +1 pour la colonne des rangs des covariables
  po <- SF$poids_attendus
  if ((length(v)+1) != (ncol(po))) browser()
  poids <- matrix(NA,nrow=nrow(po),ncol=AS$nv+1)
  poids[,1] <- po[,1]
  poids[,(v)] <- po[,-1]
  AS$poids_attendus <- poids
  return(AS)
}