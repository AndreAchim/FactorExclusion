factor_exclusion <- function(AS,AB=NULL){
  # AB contient les rangs des deux covariables dans AS$GS
  # ajoute à AS les champs $excl_AB donnant les variables représentant le facteur exclu
  # $excl_cont les contrastes de chacune des variables hors du facteur de A et B
  # $excl_crit (qui est la val. abs de la corrrel du contraste avec le témoin), $excl_poids
  # $excl_prob, $excl_cibles
  AS$excl_AB <- as.vector(AB)
  noms <- colnames(AS$GS)
  nn <- nchar(noms)
  mn <- nn[length(nn)]
  if (mn==1)
  AS$excl_AB <- as.vector(AB)
  noms <- colnames(AS$GS)
  nn <- nchar(noms)
  mn <- nn[length(nn)]
  if (mn==1)
    AS$excl_cibles <- setdiff(AS$pertinent,AB)
  else
    AS$excl_cibles <- which(nn==mn)
  cibles <- AS$excl_cibles
  ncibl <- str_split(noms,"")
  for (k in 1:2){
    nm2 <- noms[AB[k]]
    pos <- str_locate(noms[cibles],nm2)[,1]
    ciblesOK <- cibles[is.na(pos)]
    out <- exclude_fact(AS,AB[k],AB[-k],ciblesOK)
    AS$excl_poids[[k]] <- out$poids
    AS$excl_crit[[k]] <- out$crit 
    AS$excl_prob[[k]] <- out$prob 
    AS$excl_corr[[k]] <- out$corr
    nm <- paste0(noms[ciblesOK],nm2)
    colnames(out$contrasts) <- nm
    AS$GS <- cbind(AS$GS,out$contrasts)
  }
  return(AS)
}

exclude_fact <- function(AS,ote,temoin,cibles){
  # procède à l'annulation de toutes les variables de AS$GS autres que ote et temoin
  # retourne une liste incluant le contraste et le critère
  nc <- nrow(AS$GS)  # longueur des scores de chaque contraste
  nv <- length(cibles)
  con <- matrix(0,nc,nv)
  poids <- numeric(nv)  # cela laissera des 0 pour les variables pas dans AS$pertinent
  prob <- numeric(nv)
  corCT <- numeric(nv)  # correl du contraste avec le témoin
  # Saturations initiales estimées à partir de la paire (ote, temoin) — une seule fois
  sat_init <- sat_sur_facteur(AS$R, AS$N, c(ote, temoin))
  for (v in 1:nv){
    DatCol <- AS$GS[,c(cibles[v],ote,temoin)]
    # po_init = x/a_hat ; repli sur NULL si estimation invalide
    # x_v <- sat_init$sat[cibles[v]]
    x_v <- sat_init$sat[v]
    po_init <- NULL # if (!is.na(x_v) && sat_init$a_hat != 0) x_v / sat_init$a_hat else NULL
    out <- optim_exclusion(DatCol, po_init = po_init)
    corCT[v] <- out$crit
    prob[v] <- probCorr(out$crit,AS$N)$p
    poids[v] <- out$poids
    con[,v] <- out$contrast
  }
  return(list(poids=poids,corCT=corCT,prob=prob,contrasts=con))
}

optim_exclusion <- function(DatCol, po_init = NULL){
  # DatCol[,1:3] tiré de  AS$GS[,c(cibles[v],ote,temoin)]
  # col 1 est la variable à traiter
  # col 2 est la covariable pour exclure le facteur
  # col 3 est la variable témoin
  # le contraste est -v1+p*v2 où p minimise sa corrélation absolue avec v3
  # po_init : valeur initiale fournie par sat_sur_facteur (si NULL, calcul local)
  if (is.null(po_init)) {
    Rv <- t(DatCol[,2:3]) %*% DatCol[,1]
    Rab <- (t(DatCol[,2]) %*% DatCol[,3])[1]
    if (prod(c(Rv,Rab))<0)
      po <- 0
    else {
      Rva <- Rv[1]
      Rvb <- Rv[2]
      vv <- sqrt(Rva * Rvb / Rab) * sign(Rva)
      aa <- sqrt(Rab * Rva / Rvb)    # aa toujours positif
      if (is.na(aa) || is.na(vv)) browser()
      po <- vv / aa
    }
  } else {
    po <- po_init
  }
  out <- optim(po,crit_abs_cor,gr=NULL,DatCol,method="BFGS")
  contrast <- DatCol[,1:2] %*% c(-1,out$par)
  return(list(crit=out$value,poids=out$par,contrast=contrast))
}

crit_abs_cor <- function(po,Dat) {
  contrast <- sc1(Dat[,1:2] %*% c(-1,po))
  cor <- t(contrast) %*% Dat[,3]
  return(abs(cor))
}
