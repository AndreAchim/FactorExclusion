test_proj <- function(AS,cx,cy,srce=NULL){
  # cx et cy sont des rangs de variables; v_excl donne la paire d'exclusion de dimension. 
  #>>>> Paire devrait devenir un vecteur qui peut contenir plus que 2 variables exclusives au facteur à exclure ###
  # Forme pour chaque cas les contrastes AS$excl_weig[v_excl[1]][c(cx,cy)] et AS$excl_weig[v_excl[2]][c(cx,cy)]
  # Moyenne pour chaque cas les produits (cx,cy) de listes différentes
  # if (is.null(dim(AS$zdat)))  #
  #   AS <- cree_zdat(AS)       #
  browser()
  v_excl <- AS$excl_var
  ne <- length(v_excl)
  x <- which(cx==AS$excl_cibles)
  y <- which(cy==AS$excl_cibles)
  if (is.null(x) || is.null(y)) error("cx et cy doivent être des rangs de variables inclus dans AS$excl_cibles")
  XX <- matrix(NA,nrow=AS$nv,ncol=ne+1)
  YY <- matrix(NA,nrow=AS$nv,ncol=ne+1)
  XX[,1] <- AS$GS[,cx]
  YY[,1] <- AS$GS[,cy]
  browser()
  for (j in 1:ne){
    XX[,j+1] <- AS$GS[,x]  # pas génial de mettre x aux deux dernières positions`
    YY[,j+1] <- AS$GS[,y]
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
        if (r[s]>=1) browser()
        tr <- r[s]*sqrt(dl/(1-r[s]*r[s]))
        AA <- warnings();if(length(AA)>0) browser()
        pr[s] <- 2*pt(-abs(tr),dl)
      }
  }
  # if (any(sign(CR)!=sign(CR[1]))){
  #   pr <- rep(.5,ne*ne1)
  #   P <- .5
  # } else
  mr <- mean(atanh(r))
  # t <- mr*sqrt(dl/(1-mr*mr))
  # P <- 2*pt(-abs(t),dl)                  # CorrAbs
#  P <- 2*pnorm(-abs(mean(qnorm(pr))))
  P <- 2*pnorm(-abs(mr*sqrt((dl-1)))) # ProbTrans
  #  return(c(CR,r,pr,AS$satur_var))
  # C1x <-  AS$zdat[,cx] - AS$excl_weig[[1]][x]*AS$zdat[,AS$excl_var[1]]
  # C2y <-  AS$zdat[,cy] - AS$excl_weig[[2]][y]*AS$zdat[,AS$excl_var[2]]
  # sc <- C1x * AS$zdat[,cy] + C2y * AS$zdat[,cx]
  # t <- t.test(sc)
  # P <- list(P,t$p.value,t$estimate)
  # if (cx==7 && cy==8 && t$p.value<.05) browser() #print(t$p.value)
  return(list(CR=CR,r=r,pr=pr,C=mean(CR),R=tanh(mr),P=P))
}
