factor_exclusion <- function(AS,var=NULL){
  # var est un rang de variable dans AS$dat qui doit appartenir à une paire auto-annulante
  # si absent, la première variable de AS$Z est utilisée; c'est une variable avec la meilleure annulation
  # ajoute à AS les champs excl_AB donnant les deux variables représentant le facteur exclu
  # A_excl et B_excl les contrastes de chacune des variables hors du facteur de A et B
  # lorsque la variable est partiellement annulée par A ou par B, l'autr servant de témoin de l'annulation
  if (is.null(var)) var <- AS$Z[1,1]
  # trouver la grappe qui contient var
  gr <- AS$VG[[1]]$Gr
  AB <- NULL
  for (k in 1:length(gr)){
    if (any(gr[[k]]==var)){
      AB <- gr[[k]]
      break 
    }
  }
  if (is.null(AB)){
    browser()
    stop("La variable désignée n'a pas d'annulation de son signal")
    }
  AS$excl_AB <- AB
  AS$util <- setdiff(AS$pertinent,AB)
  n <- AS$N-1
#  browser()
  for (k in 1:length(AB)){
    out <- exclude_fact(AS$GS[,AS$util],AS$GS[,AB[k]],AS$GS[,AB[-k]])
    AS$excl_we[[k]] <- out$weights
    AS$excl_cr[[k]] <- out$crit * n 
    AS$excl_co[[k]] <- out$contrasts 
    AS$excl_pr[[k]] <- out$proj
#    AS$excl[[k]] <- out
  }
  return(AS)
}

exclude_fact <- function(CIBLES,OTE,TEMOINS){
# procède à l'annulation de toutes les variables de CIBLES
# avec la variable OTE, avec les variables de TEMOINS comme temoins
# retourne une liste incluant le contraste et le  (pas multiplié par (N-1))
  if (is.matrix(CIBLES)){
    nv <- ncol(CIBLES)
    nc <- nrow(CIBLES)
  } else {
    nv <- 1
    nc <- length(CIBLES)
  }
  if (is.matrix(TEMOINS))
    nt <- ncol(TEMOINS)
  else
    nt <- 1
  con <- matrix(0,nc,nv)
  wei <- numeric(nv)  # cela laissera des 0 pour les variables pas dans AS$pertinent
  cri <- numeric(nv)
  pro <- matrix(0,nt,nv)
  for (v in 1:nv){
    CiblOte <- cbind(CIBLES[,v],OTE)
    p <- optimize(pair_cancel,c(-1000,1000),CiblOte,TEMOINS)
    wei[v] <- p$minimum
    out <- pair_cancel(wei[v],cbind(CIBLES[,v,drop = FALSE],OTE),TEMOINS,as_list=TRUE)
    con[,v] <- out$contrast
    cri[v] <- out$crit
    pro[,v] <- out$proj
  }
  return(list(weights=wei,crit=cri,contrasts=con,proj=pro))
}

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