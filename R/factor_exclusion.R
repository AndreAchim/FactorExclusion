factor_exclusion <- function(AS,var=NULL){
  # var est un vecteur de rangs de variables dans AS$GS toutes exclusives à un même facteur
  # si absent, la paire qui s'annule le mieux est utilisée
  # ajoute à AS les champs excl_AB donnant les variables représentant le facteur exclu
  # A_excl et B_excl les contrastes de chacune des variables hors du facteur de A et B
  # lorsque la variable est partiellement annulée par A ou par B, l'autre servant de témoin de l'annulation
  # AS$poids_attendus, si présent, sert à remplacer l'optimisation
  if (is.null(var)){
    p <- which.max(AS$Prob)
    AB <- AS$Cpaire[,p]
  } else if (length(AB)==1) {
    p <- c(which(AS$Cpaires[1,]==AB),which(AS$Cpaires[2,]==AB))
    q <- which(AS$Prob[p]>.05)
    if (length(q)==0)
      stop("Toutes les annulation de la variable désignée ont p<.05")
    AB <- unique(as.vector(AS$Cpaires[,p[q]]))
  } else
    AB <- var
  AS$excl_AB <- as.vector(AB)
  AS$util <- AS$pertinent[-AB]
  po <- NA
  poids <- AS$poids_attendus   # NULL si le champ n'existe pas
  for (k in 1:length(AB)){
    if (!is.null(poids)){
      j <- which(poids[,1]==AB[k])
      po <- poids[j,-c(1,2,3)]
    }
    out <- exclude_fact(AS,AB[k],AB[-k],AS$util,poids=po)
#    out <- exclude_fact(AS$GS[,AS$util],AS$GS[,AB[k]],AS$GS[,AB[-k]],poids=po)
    AS$excl_weig[[k]] <- out$weights
    AS$excl_crit[[k]] <- out$crit 
    AS$excl_cont[[k]] <- out$contrasts 
    AS$excl_corr[[k]] <- out$corr
  }
#  browser()
  return(AS)
}

exclude_fact <- function(AS,ote,temoin,cibles,poids=NA){
# procède à l'annulation de toutes les variables de AS$GS autres que ote et temoin
# retourne une liste incluant le contraste et le critère
  nc <- nrow(AS$GS)  # longueur des scores de chaque contraste
  nv <- length(as.vector(cibles))
  nt <- length(as.vector(temoin))
  con <- matrix(0,nc,nv)
  wei <- numeric(nv)  # cela laissera des 0 pour les variables pas dans AS$pertinent
  cri <- numeric(nv)
#  pro <- matrix(0,nt,nv)
  cor <- matrix(0,nt,nv)
  for (v in 1:nv){
    DatCol <- AS$GS[,c(cibles[v],ote,temoin)]
    if (is.na(poids[v])){
      out <- crit_paire_initiale(AS,1,2,DatCol)
      #      out <- annule_paire(DatCol,1,2)
      # p <- optimize(pair_cancel,c(-1000,1000),CiblOte,TEMOINS)
      wei[v] <- out$po
#      browser()
    }
    else
      wei[v] <- poids[v]
    if (is.na(wei[v])) {
      out <- pair_cancel(wei[v],cbind(CIBLES[,v,drop = FALSE],OTE),TEMOINS,as_list=TRUE)
#browser()
    }
    con[,v] <- out$contrast
    cri[v] <- out$crit
#    pro[,v] <- out$proj[out$proj!=0]
    cor[,v] <- out$corr
  }
#  return(list(weights=wei,crit=cri,contrasts=con,proj=pro,corr=cor))
  return(list(weights=wei,crit=cri,contrasts=con,corr=cor))
}


# exclude__fact <- function(CIBLES,OTE,TEMOINS,poids=NA){
#   # procède à l'annulation de toutes les variables de CIBLES
#   # avec la variable OTE, avec les variables de TEMOINS comme temoins
#   # retourne une liste incluant le contraste et le critère
#   if (is.matrix(CIBLES)){
#     nv <- ncol(CIBLES)
#     nc <- nrow(CIBLES)
#   } else {
#     nv <- 1
#     nc <- length(CIBLES)
#   }
#   if (is.matrix(TEMOINS))
#     nt <- ncol(TEMOINS)
#   else
#     nt <- 1
#   con <- matrix(0,nc,nv)
#   wei <- numeric(nv)  # cela laissera des 0 pour les variables pas dans AS$pertinent
#   cri <- numeric(nv)
# #  pro <- matrix(0,nt,nv)
#   cor <- matrix(0,nt,nv)
#   for (v in 1:nv){
#     DatCol <- cbind(CIBLES[,v],OTE,TEMOINS)
#     browser()
#     if (is.na(poids[v])){
#       out <- crit_paire_initiale(AS,1,2,aaa)
# #      out <- annule_paire(DatCol,1,2)
#       # p <- optimize(pair_cancel,c(-1000,1000),CiblOte,TEMOINS)
#       wei[v] <- out$po
#       browser()
#     }
#     else
#       wei[v] <- poids[v]
#     if (is.na(wei[v])) {
#       out <- pair_cancel(wei[v],cbind(CIBLES[,v,drop = FALSE],OTE),TEMOINS,as_list=TRUE)
#     }
#     con[,v] <- out$contrast
#     cri[v] <- out$crit
# #    pro[,v] <- out$proj[out$proj!=0]
#     cor[,v] <- out$corr
#   }
# #  return(list(weights=wei,crit=cri,contrasts=con,proj=pro,corr=cor))
#   return(list(weights=wei,crit=cri,contrasts=con,corr=cor))
# }

corr_excl <- function(AS){
  ote <- AS$excl_AB
  var <- setdiff(AS$pertinent,ote)
  k <- length(var)
  corr_brut <- matrix(0,k,k)
  corr_mitig <- corr_brut
  for (j in 1:(k-1)){
    jj <- var[j]
    for (i in (j+1):k){
      ii <- var[i]
      if (any(ote==ii) | any(ote==jj))
        corr_brut[j,i] <- NA_real_
      else
        corr_brut[j,i] <- t(AS$A_excl$contrasts[,ii]) %*% AS$B_excl$contrasts[,jj]
    }
  }
  AS$cor_exc <- list(var=var,corr_brut=corr_brut,corr_mitig=corr_mitig)
  return(AS)
}


# exclude_fact <- function(CIBLES,OTE,TEMOINS,poids=NA){
#   # procède à l'annulation de toutes les variables de CIBLES
#   # avec la variable OTE, avec les variables de TEMOINS comme temoins
#   # retourne une liste incluant le contraste et le critère
#   if (is.matrix(CIBLES)){
#     nv <- ncol(CIBLES)
#     nc <- nrow(CIBLES)
#   } else {
#     nv <- 1
#     nc <- length(CIBLES)
#   }
#   if (is.matrix(TEMOINS))
#     nt <- ncol(TEMOINS)
#   else
#     nt <- 1
#   con <- matrix(0,nc,nv)
#   wei <- numeric(nv)  # cela laissera des 0 pour les variables pas dans AS$pertinent
#   cri <- numeric(nv)
#   pro <- matrix(0,nt,nv)
#   cor <- matrix(0,nt,nv)
#   for (v in 1:nv){
#     CiblOte <- cbind(CIBLES[,v],OTE)
#     if (is.na(poids[v])){
#       browser()
#       p <- optimize(pair_cancel,c(-1000,1000),CiblOte,TEMOINS)
#       wei[v] <- p$minimum
#     }
#     else
#       wei[v] <- poids[v]
#     if (is.na(wei[v]))
#       out <- pair_cancel(wei[v],cbind(CIBLES[,v,drop = FALSE],OTE),TEMOINS,as_list=TRUE)
#     con[,v] <- out$contrast
#     cri[v] <- out$crit
#     pro[,v] <- out$proj
#     cor[,v] <- out$corr
#   }
#   return(list(weights=wei,crit=cri,contrasts=con,proj=pro,corr=cor))
# }

