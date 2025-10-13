crit_paire_initiale <- function(AS,a,b,DT=NULL){
  # DT, si présent est soit un remplacement pour AS$Gs, 
  # soit la liste des rangs de colonnes de AS$GS à utiliser comme témoins
  inv <- FALSE  # être sur que c'est défini
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
    tem <- tem[-c(a,b)]
    TEMOINS <- as.matrix(AS$GS[,tem])   # pour quand ce serait une seule colonne
  }
  #if (a==5 && b==8) browser()
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
    cr <- f_crit(-p,col,TEMOINS)
    if (cr < p_c[3,2]) {
      p_c[3,] <- c(-p,cr)
      p_c <- p_c[order(p_c[,2]),]
    }
    cr <- f_crit(p,col,TEMOINS)
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
    out <- optimize(f_crit,limites,col,TEMOINS)
  } else {
    pas <- (limites[2]-limites[1])/20
    li <- seq(limites[1],limites[2],pas)
    cr <- li
    for (p in 1:length(li))
      cr[p] <- f_crit(li[p],col,TEMOINS)
    limites <- limites_balayage(li,cr)
    #if (!is.null(DT)) browser()
    out <- optimize(f_crit,limites,col,TEMOINS)
  }
  po <- out$minimum
  crit <- out$objective
  pr <- prodCorr(crit,AS$N)
  prob <- 1-(1-pr$p)^length(tem)
  cc <- 1:2
  if (inv){
    po <- 1/po
    cc <- c(2,1)
  }
  corr <- t(sc1(col[,cc] %*% c(1,-po))) %*%  TEMOINS
  contrats <- col %*% c(1,-po)
  return(list(crit=crit,po=po,prob=prob,corr=corr,contrast=contrats))
}


f_crit <- function(p,col,TEMOINS)
  max(abs(t(sc1(col %*% c(1,-p))) %*%  TEMOINS))

limites_balayage <- function(li,cr){
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

# optise <- function(fun,lim,A,B){
#   if (!is.numeric(lim) || length(lim)!= 2 || lim[2]<=lim[1]) browser
#   print(lim)
#   out <- optimize(fun,lim,A,B) 
# }

# annule_paire <- function(dat,a,b,N=NULL){
#   # dat est (N,nv) ou mieux, (nv,nv) obtenu par chol(corr(data))
#   # Les colonnes de dat doivent avoir une somme de carrés de 1.0.
#   # Produit le meilleur sc1(contraste dat[,c(a,b)] %*% c(1,-p)) pour
#   # minimiser max(abs(proj)) appliqué aux projections des variables 
#   # restantes sur ce contraste normalisé.
#   # ?? Retourne le critère, sa probabilité, le poids, les projections et le contraste optimal.
#   # le critère n'est-il pas -log10(prob) ??
#   nv <- ncol(dat)
#   ns <- nrow(dat)
#   ab <- cbind(dat[,a],dat[,b])  # c'est b qui est pondérée
#   temoins <- dat[,-c(a,b)]
#   nt <- ncol(as.matrix(temoins))
#   
#   opt <- crit_paire_initiale(AS,a,b)
#   
#   limites <- prep_limites_optimize(ab,temoins,signe)
#   # co <- t(dat[,a]) %*% dat[,b]  # corrélation des deux variables
#   # out <- prodCorr(co,AS$N)
#   # if (out$p > .4) {  # son signe pourrait ne pas être correct: pas de logarithme
#   #   signe <- 0
#   #   lp <- seq(-200,200,20)
#   #   browser()
#   # }
#   # else {
#   #   signe <- sign(co)
#   #   lp <- seq(-3,3,.3)
#   #   elp <- exp(lp)
#   # }
#   # # d'abord balayer la gamme raisonnable du poids recherché
#   # np <- length(lp)
#   # cr <- rep(0,np)  # critère pour chaque lp
#   # for (l in 1:np)
#   #   cr[l] <- crit_annule_paire(lp[l],ab,temoins,signe)
#   # loc <- minima_locaux(cr)
#   # po <- loc$rang
#   # lpo <- length(po)
#   # if (signe==0){
#   #   if (lpo == 0)  # pas de minimum dans la gamme balayée
#   #     if (cr[1] < cr[np])
#   #       limites <- c(-5000,-200)
#   #   else
#   #     limites <- c(200,5000)
#   # }
#   # else 
#   #   if (lpo == 0){  # pas de minimum dans la gamme balayée
#   #     if (cr[1] < cr[np])
#   #       limites <- c(-5,-2.7)
#   #     else
#   #       limites <- c(2.7,5)
#   #   }
#   # else {
#   #   if (lpo > 1) # plusieurs minima, il faut choisir
#   #     po <- po[which.min(loc$val)]
#   #   pog <- po-1  # gérer d'éventuelles égalités au minimum
#   #   while (pog>1 && cr[pog]==cr[po]) pog <- pog-1
#   #   pod <- po+1
#   #   while (pod<np && cr[pod]==cr[po]) pod <- pod+1
#   #   limites <- lp[c(pog,pod)]
#   #   if (pog >= pod) browser()
#   # }
#   out <- optimize(crit_annule_paire,limites,ab,temoins,signe)
#   resultat <- crit_annule_paire(out$minimum,ab,temoins,signe,as_list=TRUE)
#   proj <- rep(0,nv)
#   if (is.null(resultat$proj)) browser()
#   proj[-c(a,b)] <- resultat$proj
#   #    browser()
#   resultat$proj <- proj
#   resultat$po <- signe * exp(out$minimum)
# }
# if (!is.null(N)){
#   r <- resultat$crit
#   t <- r/sqrt((1-r*r)/(N-2))
#   p <- 2*pt(-abs(t),N-2)
#   if (nt>1)
#     p <- max(1-(1-p)^nt,1e-15)
#   resultat$prob <- p
#   return(resultat)
# }


#   lmt <- matrix(c(-5000,-175,175,5000,-5,-2.7,2.7,5),nrow=2)
#   co <- t(ab[,1]) %*% ab[,2]  # corrélation des deux variables
#   out <- prodCorr(co,AS$N)
#   if (out$p > .4) {  # son signe pourrait ne pas être correct: pas de logarithme
#     signe <- 0
#     lp <- seq(-200,200,20)
#     browser()
#   }
#   else {
#     signe <- sign(co)
#     lp <- seq(-3,3,.3)
#     #  elp <- exp(lp)
#   }
#   # d'abord balayer la gamme raisonnable du poids recherché
#   np <- length(lp)
#   cr <- rep(0,np)  # critère pour chaque lp
#   for (l in 1:np)
#     cr[l] <- crit_annule_paire(lp[l],ab,temoins,signe)
#   loc <- minima_locaux(cr)
#   po <- loc$rang
#   lpo <- length(po)
#   if (lpo == 0) { # pas de minimum dans la gamme balayée
#     if (cr[1] < cr[np])
#       if (signe==0){
#         limites <- c(-5000,-200)
#         else
#           limites <- c(-5,-2.7)
#       }
#     else
#       limites <- c(200,5000)
#   }
#   else 
#     if (lpo == 0){  # pas de minimum dans la gamme balayée
#       if (cr[1] < cr[np])
#         limites <- 
#           else
#             limites <- c(2.7,5)
#     }
#   else {
#     if (lpo > 1) # plusieurs minima, il faut choisir
#       po <- po[which.min(loc$val)]
#     pog <- po-1  # gérer d'éventuelles égalités au minimum
#     while (pog>1 && cr[pog]==cr[po]) pog <- pog-1
#     pod <- po+1
#     while (pod<np && cr[pod]==cr[po]) pod <- pod+1
#     limites <- lp[c(pog,pod)]
#     if (pog >= pod) browser()
#   }
#   out$limites <- limites
#   return(out)
# }


# crit_paire_initiale <- function(AS,a,b,DT=NULL){
#   # DT, si présent est soit un remplacement pour AS$Gs, 
#   # soit la liste des rangs de colonnes de AS$GS à utiliser comme témoins
#   if (is.matrix(DT)) {
#     col <- DT[,c(a,b)]
#     TEMOINS <- as.matrix(DT[,-c(a,b)])
#   } else if (is.numeric(DT)) {
#     col <- AS$GS[,c(a,b)]
#     TEMOINS <- AS$GS[,DT]
#   } else {
#     R <- abs(AS$R)
#     diag(R) <- 0
#     rA <- which.max(R[a,])
#     rB <- which.max(R[b,])
#     inv <- R[a,rA] < R[b,rB]
#     if (inv) {
#       col <- AS$GS[,c(b,a)]
#     } else {
#       col <- AS$GS[,c(a,b)]
#     }
#     tem <- AS$pertinent
#     tem <- tem[-c(a,b)]
#     TEMOINS <- as.matrix(AS$GS[,tem])   # pour quand ce serait une seule colonne
#   }
#   #if (a==5 && b==8) browser()
#   p <- .001
#   p_c <- matrix(0,nrow=3,ncol=2)
#   #  if (is.matrix(DT)) browser()
#   p_c[1,2] <- max(abs(t(col[,1]) %*% TEMOINS))
#   p_c[2,1] <- -p
#   p_c[2,2] <- max(abs(t(sc1(col %*% c(1,-p))) %*%  TEMOINS))
#   p_c[3,1] <- p
#   p_c[3,2] <- max(abs(t(sc1(col %*% c(1,p))) %*%  TEMOINS))
#   p_c <- p_c[order(p_c[,2]),]
#   while (p != 0 && p<30){
#     p <- 2 * p
#     cr <- f_crit(-p,col,TEMOINS)
#     if (cr < p_c[3,2]) {
#       p_c[3,] <- c(-p,cr)
#       p_c <- p_c[order(p_c[,2]),]
#     }
#     cr <- f_crit(p,col,TEMOINS)
#     if (cr < p_c[3,2]) {
#       p_c[3,] <- c(p,cr)
#       p_c <- p_c[order(p_c[,2]),]
#     }
#     if (order(p_c[,1])[2]==1)
#       p <- 0
#   }
#   if (p==0){
#     limites <- sort(p_c[2:3,1])
#   } else {
#     u <- p_c[1,1]
#     d <- 2*min(abs(p_c[,1]))
#     limites <- c(u-d,u+d)
#     # if (a==5 && b==8) browser()
#   }
#   dlim <- limites[2]-limites[1]
#   if (dlim<=.002) {
#     browser()
#     out <- optimize(f_crit,limites,col,TEMOINS)
#   } else {
#     pas <- (limites[2]-limites[1])/20
#     li <- seq(limites[1],limites[2],pas)
#     cr <- li
#     for (p in 1:length(li))
#       cr[p] <- f_crit(li[p],col,TEMOINS)
#     loc <- minima_locaux(cr)
#     po <- loc$rang
#     lpo <- length(po)
#     if (lpo > 1) # plusieurs minima, il faut choisir
#       po <- po[which.min(loc$val)]
#     if (po>0){
#       pog <- po-1  # gérer d'éventuelles égalités au minimum
#       while (pog>1 && cr[pog]==cr[po]) pog <- pog-1
#       pod <- po+1
#       while (pod<np && cr[pod]==cr[po]) pod <- pod+1
#       limites <- li[c(pog,pod)]
#     } else {
#       lcr <- length(cr)
#       if (cr[1]<cr[lcr]){
#         limites <- li[1:2]
#       } else
#         limites <- li[c(lcr-1,lcr)]
#     }
#     browser()
#     out <- optimize(f_crit,limites,col,TEMOINS)
#   }
#   po <- out$minimum
#   crit <- out$objective
#   pr <- prodCorr(crit,AS$N)
#   prob <- 1-(1-pr$p)^length(tem)
#   if (inv)
#     po <- 1/po
#   corr <- t(sc1(col %*% c(1,-p))) %*%  TEMOINS
#   contrats <- sc1(col %*% c(1,-p))
#   return(list(crit=crit,po=po,prob=prob,corr=corr,contrast=contrats))
# }


