factor_exclusion <- function(AS,var=NULL){
  # var est un rang de variable dans AS$dat qui doit appartenir à une paire auto-annulante
  # si absent, la première variable de AS$Z est utilisée; c'est une variable avec la meilleure annulation
  # ajoute à AS les champs excl_AB donnant les deux variables représentant le facteur exclu
  # A_excl et B_excl les contrastes de chacune des variables hors du facteur de A et B
  # lorsque la variable est partiellement annulée par A ou par B, l'autr servant de témoin de l'annulation
  if (is.null(var)) var <- AS$Z[1,1]
  fct <- AS$VG[[1]]$Fct
  factr <- which(fct[var,] != 0)[[1]]
  AB <- which(fct[,factr] != 0)
  AS$excl_AB <- AB
  AS$A_excl <- exclude_fact(AS,AB[1],AB[2])
  AS$B_excl <- exclude_fact(AS,AB[2],AB[1])
  return(AS)
}

exclude_fact <- function(AS,A,B){
  con <- matrix(0,AS$nv,AS$nv)
  wei <- numeric(AS$nv)  # cela laissera de 0 pour les variables pas dans AS$pertinent
  cri <- numeric(AS$nv)
  pro <- numeric(AS$nv)
  n=AS$N-1
  for (v in AS$pertinent){
      p <- optimize(pair_cancel,c(-100,100),AS$GS[,c(A,v)],AS$GS[,B])
      wei[v] <- p$minimum
      out <- pair_cancel(wei[v],AS$GS[,c(A,v)],AS$GS[,B],as_list=TRUE)
      con[,v] <- out$contrast
      cri[v] <- out$crit * n
      pro[v] <- out$proj
  }
  return(list(weights=wei,crit=cri,contrasts=con,proj=pro))
}

corr_excl <- function(AS){
  var <- setdiff(AS$pertinent,AS$excl_AB)
  k <- length(var)
  corr_brut <- matrix(0,k,k)
  corr_mitig <- corr_brut
  for (j in 1:(k-1)){
    jj <- var[j]
    for (i in (j+1):k){
      ii <- var[i]
      corr_brut[j,i] <- t(AS$A_excl$contrasts[,ii]) %*% AS$B_excl$contrasts[,jj]
    }
  }
  AS$cor_exc <- list(var=var,corr_brut=corr_brut,corr_mitig=corr_mitig)
  return(AS)
}