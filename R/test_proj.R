test_proj <- function(AS,cx,cy,srce=NULL){
  # cx et cy sont des rangs de variables; v_excl donne la paire d'exclusion de dimension. 
  #>>>> Paire devrait devenir un vecteur qui peut contenir plus que 2 variables exclusives au facteur à exclure ###
  # Forme pour chaque cas les contrastes AS$excl_weig[v_excl[1]][c(cx,cy)] et AS$excl_weig[v_excl[2]][c(cx,cy)]
  # Moyenne pour chaque cas les produits (cx,cy) de listes différentes
  ### ci-dessous, description à actualiser
  # En fait la somme et détermine son rang parmi 1000 cas où la moitié des signes seraient inversés
  # Cela est fait en comparant la somme de sous-ensembles aléatoires avec 0
  v_excl <- AS$excl_var
  ne <- length(v_excl)
  x <- which(cx==AS$excl_cibles)
  y <- which(cy==AS$excl_cibles)
  if (is.null(x) || is.null(y)) error("cx et cy doivent être des rangs de variables inclus dans AS$util")
  XX <- matrix(NA,nrow=AS$nv,ncol=ne+1)
  YY <- matrix(NA,nrow=AS$nv,ncol=ne+1)
  XX[,1] <- AS$GS[,cx]
  YY[,1] <- AS$GS[,cy]
  for (j in 1:ne){
    XX[,j+1] <- AS$excl_cont[[j]][,x]
    YY[,j+1] <- AS$excl_cont[[j]][,y]
  }
  ne1 <- ne+1
  CR <- rep(NA,ne*ne1)
  r <- rep(NA,ne*ne1)
  pr <- rep(NA,ne*ne1)
  s <- 0
  dl <- AS$N-2
  for (j in 1:ne1){
    for (k in 1:ne1)
      if (k!=j){
        s <- s+1
        CR[s] <- t(XX[,j]) %*% YY[,k]
        r[s] <- t(sc1(XX[,j])) %*% sc1(YY[,k])
        tr <- r[s]*sqrt(dl/(1-r[s]*r[s]))
        pr[s] <- 2*pt(-abs(tr),dl)
      }
  }
  # if (any(sign(CR)!=sign(CR[1]))){
  #   pr <- rep(.5,ne*ne1)
  #   P <- .5
  # } else
  mr <- mean(r)
  t <- mr*sqrt(dl/(1-mr*mr))
  P <- 2*pt(-abs(t),dl)
  # P <- 2*pnorm(-abs(mean(qnorm(pr))))
  #  return(c(CR,r,pr,AS$satur_var))
  return(list(CR=CR,r=r,pr=pr,C=mean(CR),R=mean(r),P=P))
}
