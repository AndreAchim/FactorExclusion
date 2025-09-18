test_proj <- function(AS,cx,cy,paire=c(1,2)){
  # cx et cy sont des rangs de variables; paire donne la paire d'exclusion de dimension. 
  #>>>> Paire devrait devenir un vecteur qui peut contenir plus que 2 variables exclusives au facteur à exclure ###
  # Forme pour chaque cas les contrastes AS$excl_weig[paire[1]][c(cx,cy)] et AS$excl_weig[paire[2]][c(cx,cy)]
  # Moyenne pour chaque cas les produits (cx,cy) de listes différentes
  ### ci-dessous, description à actualiser
  # En fait la somme et détermine son rang parmi 1000 cas où la moitié des signes seraient inversés
  # Cela est fait en comparant la somme de sous-ensembles aléatoires avec 0
  if (is.null(dim(AS$zdat)))
    AS <- cree_zdat(AS)
  x <- which(cx==AS$util)
  y <- which(cy==AS$util)
  if (is.null(x) || is.null(y)) error("cx et cy doivent être des rangs de variables inclus dans AS$util")
  # # à compléter
  #   CC <- 0
  #   for (un in paire)
  #     for (deux in paire)
  #       if (un!=deux){
  #         cc1 <- AS$zdat[,cx] - AS$excl_weig[[un]][x]*AS$zdat[,AS$excl_AB[un]]
  #         cc2 <- AS$zdat[,cx] - AS$excl_weig[[un]][x]*AS$zdat[,AS$excl_AB[un]]
  #       }
  
  # xy <- c(x,y)
  # for (k in paire){
  #   po <- AS$excl_weig[[k]][xy]
  #   Cont <- AS$zdat[,paire,]
  C1x <-  AS$zdat[,cx] - AS$excl_weig[[1]][x]*AS$zdat[,AS$excl_AB[1]]
  C1y <-  AS$zdat[,cy] - AS$excl_weig[[1]][y]*AS$zdat[,AS$excl_AB[1]]
  C2x <-  AS$zdat[,cx] - AS$excl_weig[[2]][x]*AS$zdat[,AS$excl_AB[2]]
  C2y <-  AS$zdat[,cy] - AS$excl_weig[[2]][y]*AS$zdat[,AS$excl_AB[2]]
  #  CC <- (C1x * C2y + C2x * C1y)/2
  N=length(C1x)
  out <- vector()
  for (fois in 1:2){
    if (length(out)){
      CC <- C2x * C1y
      r <- sum(CC) / sqrt((t(C1y) %*% C1y) * (t(C2x) %*% C2x))
    }
    else {
      CC <- C1x * C2y 
      r <- sum(CC) / sqrt((t(C1x) %*% C1x) * (t(C2y) %*% C2y))
    }
    #  proj_moy <- sum(CC)
    # approche par inversion de signes abandonnée au profit du t de Student sur une moyenne
    # proj_crit <- rep(0,1000)
    # for (k in 1:1000)
    #   proj_crit[k] <- sum(sample(CC,floor(AS$N/2)))
    # if (proj_moy>0)
    #   proj_prob <- sum(proj_crit<0)
    # else
    #   proj_prob <- sum(proj_crit>0)
    # proj_prob <- (1+proj_prob)/1001
    # return(list(proj_moy=proj_moy,proj_prob=proj_prob, proj_crit=proj_crit,Prod=CC))
    
    # r12 <- (t(C1x) %*% C2y) /sqrt((t(C1x) %*% C1x) * (t(C2y) %*% C2y))
    # r21 <- (t(C1y) %*% C2x) /sqrt((t(C1y) %*% C1y) * (t(C2x) %*% C2x))
    # r <- (r12+r21)/2
    dl <- N-3
    tr <- r*sqrt(dl/(1-r*r))
    pr <- 2*pt(-abs(tr),dl)
    tt1 <- t.test(CC)
    t1e <- tt1$estimate*N
    t1t <- tt1$statistic
    t1p <- tt1$p.value
    p1b <- 2*pt(-abs(t1t),dl)
    tt2 <- t.test(sign(CC) * sqrt(abs(CC)))
    t2e <- tt2$estimate*sqrt(N)
    t2t <- tt2$statistic
    t2p <- tt2$p.value
    p2b <- 2*pt(-abs(t2t),dl)
    out <- c(out,r,t1e,t2e,tr,t1t,t2t,pr,t1p,t2p,p1b,p2b)
  }
    return(unname(out))
}