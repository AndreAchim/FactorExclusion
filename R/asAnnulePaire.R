asAnnulePaire <- function(AS, k) {
  # AS$GS(v,v) est la transformation chol() (en triangle supérieur) de R(v,v)
  # k est un scalaire pointant dans AS$Cpaires(np,2)
  # Optimise les variables de la rangée k de AS$Cpaires pour en minimiser le signal
  # Remplit AS$Crit[k], AS$Ppaires[k], AS$Corr[,k], AS$Prob[k]
  # AS$absLogPo sert à indiquer une annulation pas réussie (poids très petit ou très grand)
  # AS$Crit[k] sera R2 de régression multple   <<<<<<<<<
  melange <- AS$Cpaires[, k]
  cible <- setdiff(AS$pertinent,melange)
  pd <- probNonDoublet(AS, melange)
  r <- AS$R[melange[1],melange[2]]  # corrélation des variables mêmes
  if (pd < 0.05) {  # si la paire a d'autres corrélations
    ap <- crit_paire_initiale(AS,melange[1],melange[2])
  } else {    # si la paire semble dépendre d'un facteur doublet
    AS$doublet <- rbind(AS$doublet, melange)  # documenter le doublet dans AS
    out <- probCorr(r,AS$N)
    cont <- AS$GS[,melange] %*% c(1,-1)
    corr <- t(sc1(cont)) %*% sc1(AS$GS[,-melange])
    R2 <- crit_R2(1,AS$GS[,melange],AS$GS[,cible])  # ???? pour variables d'un doublet ????
    prob <- prob_R2(R2,length(cible),AS$N)
    ap <- list(prob=prob,po=sign(r),crit=R2,corr=corr)
    ap$po <- sign(r)  # Pour un doublet, les 2 variables auront le même poids en val. abs.
  }
  if (abs(log(abs(ap$po)))>5){
    AS$Prob[k] <- .1^abs(log(abs(ap$po)))
  }
  else
    AS$Prob[k] <- ap$prob
  AS$Crit[k] <- ap$crit
  if (length(AS$Corr[,k]) != length(ap$corr)) browser()
  AS$Corr[,k] <- ap$corr
  if (any(is.nan(ap$corr))) browser()
  AS$Ppaires[k] <- ap$po
  if (r*ap$po>0) {
    s1 <- sqrt(r*ap$po)
    s2 <- r/s1
  } else {
    s1 <- s2 <- NA
  }
  AS$satPaires[,k] <- c(s1,s2)
  return(AS)
}
