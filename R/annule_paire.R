annule_paire <- function(dat,a,b,N=NULL){
  # dat est (N,nv) ou mieux, (nv,nv) obtenu par chol(cor(data))
  # Les colonnes de dat doivent avoir une somme de carrés de 1.0.
  # Si N est connu (rangées de dat si pas nv) ou donné dans N et
  # si la corrélation des deux variables n'est pas significative à p<.25
  # l'échec de l'annulation est rapporté sans plus de calcul.
  # Produit le meilleur sc1(contraste dat[,c(a,b)] %*% c(p,-1)) pour
  # minimiser le critère fct appliqué aux projections des variables 
  # restantes sur ce contraste normalisé.
  # Si un minimum existe, retourne le critère, le poids, les projections
  # et le contraste optimal.
  # Autrement,retourne le critère et le poids comme NaN.
  code <- 0 # code de sortie pour annulation réussie
  nv <- ncol(dat)
  ns <- nrow(dat)
  if (is.null(N) && ns>nv) N=ns
  ab <- cbind(dat[,a],dat[,b])
  temoins <- dat[,-c(a,b)]
  nt <- ncol(temoins)
  co <- t(dat[,a]) %*% dat[,b]
  if (!is.null(N)){
    z <- abs(co) * sqrt(N-1)
    if (length(abs(z)) == 0) browser()
    if (2*pnorm(-abs(z))>.01)
      return(list(code=1,crit=NaN,po=Inf,prob=NaN,proj=rep(0,nv)))
  signe <- sign(co)
  lp <- seq(-3,3,.3)
  elp <- exp(lp)
  np <- length(lp)
  cr <- rep(0,np)
  for (l in 1:np)
    cr[l] <- crit_annule_paire(lp[l],ab,temoins,signe)
  loc <- minima_locaux(cr)
  po <- loc$rang
  lpo <- length(po)
  if (lpo != 1)
    if (lpo == 0){
      po <- 0  # pas de minimum dans l'intervalle balayé
      code <- 2 + (cr[1] < cr[np])  # 2 pour min à gauche et 3 à droite
#      plot(lp,cr)
#      browser()
    }
    else    # plusieurs minima, il faut choisir
      po <- po[which.min(loc$val)]
  if (po>1 && po<np){
    pog <- po-1
    while (pog>1 && cr[pog]==cr[po]) pog <- pog-1
    pod <- po+1
    while (pod<np && cr[pod]==cr[po]) pod <- pod+1
      limites <- lp[c(pog,pod)]
    out <- optimize(crit_annule_paire,limites,ab,temoins,signe)
    resultat <- crit_annule_paire(out$minimum,ab,temoins,signe,as_list=TRUE)
    proj <- rep(0,nv)
    if (is.null(resultat$proj)) browser()
    proj[-c(a,b)] <- resultat$proj
#    browser()
    resultat$proj <- proj
    resultat$po <- exp(out$minimum)
    resultat$code <- 0
    if (!is.null(N)){
      r <- resultat$crit
      if (is.nan(r))
        resultat$prob <- NaN
      else {
        t <- r/sqrt((1-r*r)/(N-2))
        p <- 2*pt(-abs(t),N-2)
        if (nt>1)
          p <- 1-(1-p)^nt
        resultat$prob <- p
      }
    }
#    if (is.numeric(N)) resultat$z2 <- resultat$crit^2 * (N-1)
#    if (is.numeric(N)) resultat$z2 <- resultat$crit * (N-1)
    return(resultat)
  } else {
#    plot(lp,cr)
    return(list(code=code,crit=NaN,po=NaN,prob=NaN,proj=rep(0,nv)))
  }
  }
}
