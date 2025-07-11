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
  if (is.null(AB))
    stop("La variable désignée n'a pas d'annulation de son signal")
  AS$excl_AB <- AB
#  browser()
  for (k in 1:length(AB)){
    out <- exclude_fact(AS,AB[k],AB[-k])
    AS$excl_we[[k]] <- out$weigths
    AS$excl_cr[[k]] <- out$crit 
    AS$excl_co[[k]] <- out$contrasts 
    AS$excl_pr[[k]] <- out$proj
#    AS$excl[[k]] <- out
  }
  return(AS)
}

exclude_fact <- function(AS,A,B){
# procède à l'annulation de toutes les variables (sauf celles dans A ou B)
# avec la variable de rang A avec les variables de rangs B comme témoins
# retourne une liste incluant le contraste et le critère 
  con <- matrix(0,AS$nv,AS$nv)
  wei <- numeric(AS$nv)  # cela laissera de 0 pour les variables pas dans AS$pertinent
  cri <- numeric(AS$nv)
  pro <- matrix(0,length(B),AS$nv)
  n=AS$N-1
  for (v in AS$pertinent){
    if(!any(v==c(A,B))){
      p <- optimize(pair_cancel,c(-100,100),AS$GS[,c(A,v)],AS$GS[,B])
      wei[v] <- p$minimum
      out <- pair_cancel(wei[v],AS$GS[,c(A,v)],AS$GS[,B],as_list=TRUE)
      con[,v] <- out$contrast
      cri[v] <- out$crit * n
      pro[,v] <- out$proj
    }
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