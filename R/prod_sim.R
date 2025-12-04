prod_sim <- function(SF,N,nrep=100,ON=FALSE){
  # SF est une liste avec les champs $patronm $corrFct et $paires. 
  # Ce dernier champ contient une liste de paires de rangs de variables (2,k)
  # tester en dehors du facteur exprimé dans les variables 1 et 2.
  # N est la taille des échantillons et nrep le nombre de jeux de données simulées
  # En sortie, on aura une liste avec, pour chaque entrée de $paires, une matrice encore à déterminer
  if (N<0){
    N <- -N
    set.seed(N)
  }
  if (isSymmetric(SF$corrFct))
    F <- SF$patron %*% sqrtm(SF$corrFct)
  else
    F <- SF$patron %*% SF$corrFct
  # réserver l'espace
  result <- list()
  k <- 0
  for (Li in SF$paires){
    k <- k + 1
    np[k] <- ncol(Li)
    result[[k]] <- array(NA,c(np[k],20,nrep))
  }
  #  sdif <- rep(0,np[k])
  # effectuer les nrep simulations
  # options(warn=1)
  for (rp in 1:nrep){
    #   print(rp)
    out <- gen_data(F,N,ON=ON)
    dt <- out$dt
    #    if (rp==15) browser()
    srce <- ortho_norm(out$srce[,1:3])  # porter chaque colonne à une somme de carrés de 1.0
    AS <- init_SCA(dt)
    #    AS <- cree_zdat(AS)
    if (!is.null(SF$poids_attendus)) AS <- prepare_poids_fixes(AS,SF)
    #    AS <- factor_exclusion(AS,c(1,2))
    AS <- exclusions(AS,c(1,2))
    AS <- cree_zdat(AS) ####
    k <- 0
    #    a <- FALSE
    for (Li in SF$paires){
      #      a <- !a
      k <- k + 1
      for (p in 1:np[k]){
        out <- test_proj(AS,SF$paires[[k]][1,p],SF$paires[[k]][2,p],srce=srce)
        out <- c(out$CR,out$r,out$pr,AS$satur_var)
        #       if (a){
        #        aa <- test_cov_resid(AS,SF$paires[[k]][1,p],SF$paires[[k]][2,p])
        #          print(c(rp,k,p,aa,mean(out[13:18])))
        #        }
        # return(list(CR_=CR_,CR=CR,r=r,PR=PR,pr_=pr_,pr=pr,pond=pond))
        # result[[k]][p,,rp] <- c(out$CR_,mean(out$CR),mean(out$r),out$PR,out$pr_,mean(out$pr))
        #        if (length(out) != 18) browser()
        # browser()
        result[[k]][p,,rp] <- out
        # if (p==np[k] && all(sign(out[2:6])==sign(out[1]))){
        #   aa <-  exp(mean(log(out[13:18])))
        #   if (Li[1,1]==3 && aa<.05){
        #     V1 <- AS$zdat[,1]
        #     V2 <- AS$zdat[,2]
        #     V7 <- AS$zdat[,7]
        #     V8 <- AS$zdat[,8]
        #     X1 <- sc1(AS$zdat[,7]-AS$excl_w[[1]][5]*AS$zdat[,1])
        #     Y1 <- sc1(AS$zdat[,8]-AS$excl_w[[1]][6]*AS$zdat[,1])
        #     X2 <- sc1(AS$zdat[,7]-AS$excl_w[[2]][5]*AS$zdat[,2])
        #     Y2 <- sc1(AS$zdat[,8]-AS$excl_w[[2]][6]*AS$zdat[,2])
        #     AA <- t(srce[,1:3]) %*% cbind(V1,V2,V7,V8,X1,X2,Y1,Y2)
        #     #            print(sprintf("%3d %6.4f   %6.4f %6.4f %6.4f %6.4f",rp,aa,AA[1],AA[2],AA[3],AA[4]),quote=FALSE)
        #     print(sprintf("rp%d p%d Li%d_%d p%6.4f",rp,p,Li[1,p],Li[2,p],aa),quote=FALSE)
        #     print(cbind(t(AS$excl_weig[[1]][5:6]),t(AS$excl_weig[[2]][5:6]),t(AS$excl_crit[[1]][5:6]),t(AS$excl_crit[[2]][5:6])))
        #     #            if (aa<.01){
        #     print(AA)
        #     browser()
        #   }
        #    #           }
        # } else
        #   sdif[k] <- sdif[k]+1
        #        result[[k]] <- out
        #        if(rp==1) print(c(p,out))
        # }
      # browser()
      # result[[k,rp]] <- out
    }
  }
}
return(result)
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