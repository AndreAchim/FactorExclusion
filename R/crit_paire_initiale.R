crit_paire_initiale <- function(AS,a,b,DT=NULL){
  # DT, si présent est soit un remplacement pour AS$Gs, 
  # soit la liste des rangs de colonnes de AS$GS à utiliser comme témoins
  inv <- FALSE  # être sûr que c'est défini
  if (is.matrix(DT)) {
    col <- DT[,c(a,b)]
    TEMOINS <- as.matrix(DT[,-c(a,b)])
  } else if (is.numeric(DT)) {
    col <- AS$GS[,c(a,b)]
    TEMOINS <- AS$GS[,DT]
  } else {
    R <- abs(AS$R)
    diag(R) <- 0
    rA <- which.max(R[a,])
    rB <- which.max(R[b,])
    inv <- R[a,rA] < R[b,rB]
    if (inv) {
      col <- AS$GS[,c(b,a)]
    } else {
      col <- AS$GS[,c(a,b)]
    }
    tem <- AS$pertinent
    tem <- tem[-c(which(tem==a),which(tem==b))]
    TEMOINS <- as.matrix(AS$GS[,tem])   # pour quand ce serait une seule colonne
  }
  Rinv <- solve(t(TEMOINS) %*% TEMOINS)
  p <- .001
  p_c <- matrix(0,nrow=3,ncol=2)
  #  if (is.matrix(DT)) browser()
  p_c[1,2] <- max(abs(t(col[,1]) %*% TEMOINS))
  p_c[2,1] <- -p
  p_c[2,2] <- max(abs(t(sc1(col %*% c(1,-p))) %*%  TEMOINS))
  p_c[3,1] <- p
  p_c[3,2] <- max(abs(t(sc1(col %*% c(1,p))) %*%  TEMOINS))
  p_c <- p_c[order(p_c[,2]),]
  while (p != 0 && p<30){
    p <- 2 * p
    cr <- crit_R2(-p,col,TEMOINS,Rinv)
    if (cr < p_c[3,2]) {
      p_c[3,] <- c(-p,cr)
      p_c <- p_c[order(p_c[,2]),]
    }
    cr <- crit_R2(p,col,TEMOINS,Rinv)
    if (cr < p_c[3,2]) {
      p_c[3,] <- c(p,cr)
      p_c <- p_c[order(p_c[,2]),]
    }
    if (order(p_c[,1])[2]==1)
      p <- 0
  }
  if (p==0){
    limites <- sort(p_c[2:3,1])
  } else {
    pc <- p_c[,1]
    u <- pc[1]
    if (any(pc==0)) pc <- pc[-which(pc==0)]
    d <- 2*min(abs(pc))
    limites <- c(u-d,u+d)
  }
  dlim <- limites[2]-limites[1]
  if (dlim<=.002) {
    out <- optimize(crit_R2,limites,col,TEMOINS,Rinv)
  } else {
    pas <- (limites[2]-limites[1])/20
    li <- seq(limites[1],limites[2],pas)
    cr <- li  # juste pour créer cr
    for (p in 1:length(li))
      cr[p] <- crit_R2(li[p],col,TEMOINS,Rinv)
    limites <- limites_balayage(li,cr)
    out <- optimize(crit_R2,limites,col,TEMOINS,Rinv)
    # R2 <- out$objective
    # probR2 <- prob_R2(R2,ncol(TEMOINS),AS$N)
    # for (p in 1:length(li))
    #   cr[p] <- f_crit(li[p],col,TEMOINS)
    # browser()
    # limites <- limites_balayage(li,cr)
    # out <- optimize(f_crit,limites,col,TEMOINS)
    # probMaxAbs <- prodCorr(out$objective,AS$N)
  }
  po <- out$minimum
  R2 <- out$objective
  prob <- prob_R2(R2,ncol(TEMOINS),AS$N)
  # pr <- prodCorr(crit,AS$N)   % <<<<<<<<   Pas pour R2 <<<<<<<
  # prob <- 1-(1-pr$p)^length(tem)  %   <<<<<<< pas à ajuster pour R2
  cc <- 1:2
  if (inv){
    po <- 1/po
    cc <- c(2,1)
  }
  corr <- t(sc1(col[,cc] %*% c(1,-po))) %*%  TEMOINS
  contrats <- col %*% c(1,-po)
  return(list(crit=R2,po=po,prob=prob,corr=corr,contrast=contrats))
}


f_crit <- function(p,col,TEMOINS)
  max(abs(t(sc1(col %*% c(1,-p))) %*%  TEMOINS))

limites_balayage <- function(li,cr){
  np <- length(li)
  loc <- minima_locaux(cr)
  po <- loc$rang
  lpo <- length(po)
  if (lpo > 1) # plusieurs minima, il faut choisir
    po <- po[which.min(loc$val)]
  if (po>0){
    pog <- po-1  # gérer d'éventuelles égalités au minimum
    while (pog>1 && cr[pog]==cr[po]) pog <- pog-1
    pod <- po+1
    while (pod<np && cr[pod]==cr[po]) pod <- pod+1
    limites <- li[c(pog,pod)]
  } else {
    lcr <- length(cr)
    if (cr[1]<cr[lcr]){
      limites <- li[1:2]
    } else
      limites <- li[c(lcr-1,lcr)]
  }
  return(limites)
} 
