annule_paire <- function(dat,a,b,N=NULL){
  # dat est (N,nv) ou mieux, (nv,nv) obtenu par chol(cor(data))
  # Les colonnes de dat doivent avoir une somme de carrés de 1.0.
  # Produit le meilleur sc1(contraste dat[,c(a,b)] %*% c(p,-1)) pour
  # minimiser max(abs(proj)) appliqué aux projections des variables 
  # restantes sur ce contraste normalisé.
  # Retourne le critère, sa probabilité, le poids, les projections et le contraste optimal.
  nv <- ncol(dat)
  ns <- nrow(dat)
  ab <- cbind(dat[,a],dat[,b])
  temoins <- dat[,-c(a,b)]
  nt <- ncol(temoins)
  co <- t(dat[,a]) %*% dat[,b]  # corrélation des deux variables
  signe <- sign(co)
  # d'abord balayer la gamme raisonnable du log du poids recherché
  lp <- seq(-3,3,.3)
  elp <- exp(lp)
  np <- length(lp)
  cr <- rep(0,np)  # critère pour chaque lp
  for (l in 1:np)
    cr[l] <- crit_annule_paire(lp[l],ab,temoins,signe)
  loc <- minima_locaux(cr)
  po <- loc$rang
  lpo <- length(po)
  if (lpo == 0)  # pas de minimum dans la gamme balayée
    if (cr[1] < cr[np])
      limites <- c(-5,-2.7)
    else
      limites <- c(2.7,5)
  else {
    if (lpo > 1) # plusieurs minima, il faut choisir
      po <- po[which.min(loc$val)]
    pog <- po-1  # gérer d'éventuelles égalités au minimum
    while (pog>1 && cr[pog]==cr[po]) pog <- pog-1
    pod <- po+1
    while (pod<np && cr[pod]==cr[po]) pod <- pod+1
      limites <- lp[c(pog,pod)]
  }
  out <- optimize(crit_annule_paire,limites,ab,temoins,signe)
  resultat <- crit_annule_paire(out$minimum,ab,temoins,signe,as_list=TRUE)
  proj <- rep(0,nv)
  if (is.null(resultat$proj)) browser()
  proj[-c(a,b)] <- resultat$proj
#    browser()
  resultat$proj <- proj
  resultat$po <- exp(out$minimum)
  if (!is.null(N)){
    r <- resultat$crit
    t <- r/sqrt((1-r*r)/(N-2))
    p <- 2*pt(-abs(t),N-2)
    if (nt>1)
      p <- max(1-(1-p)^nt,1e-15)
    resultat$prob <- p
  }
return(resultat)
}
