asAnnulePaire <- function(AS, k) {
  # AS$GS(v,v) est la transformation chol() (en triangle supérieur) de R(v,v)
  # k est un scalaire pointant dans AS$Cpaires(np,2)
  # Optimise les variables de la rangée k de AS$Cpaires pour en minimiser le signal
  # Remplit AS$Crit[k], AS$Ppaires[k], AS$Corr[,k], AS$Prob[k]
  # AS$absLogPo sert à indiquer une annulation pas réussie (poids très petit ou très grand)
  # AS$Crit[k] sera max(abs(corrélations du contraste))
  # AS$Prob pourra être transformé copr devenir des distances pour l'analyse en grappes
  cible <- AS$pertinent
  melange <- AS$Cpaires[, k]
  pd <- probNonDoublet(AS, melange)
  if (pd < 0.05) {  # si la paire a d'autres corrélations
    ap <- crit_paire_initiale(AS,melange[1],melange[2])
  } else {    # si la paire semble dépendre d'un facteur doublet
    AS$doublet <- rbind(AS$doublet, melange)  # documenter le doublet dans AS
    r <- AS$R[melange[1],melange[2]]  # corrélation des variables mêmes
    out <- prodCorr(r,AS$N)
    cont <- AS$GS[,melange] %*% c(1,-1)
    corr <- t(sc1(cont)) %*% AS$GS[,-melange]
    ap <- list(prob=out$p,po=sign(r),crit=max(abs(corr)),corr=corr)
    ap$po <- sign(r)  # Pour un doublet, les 2 variables auront le même poids en val. abs.
    }
  if (abs(log(abs(ap$po)))>5){
    AS$Prob[k] <- .1^abs(log(abs(ap$po)))
  }
  else
    AS$Prob[k] <- ap$prob
  AS$Crit[k] <- ap$crit
  AS$Corr[,k] <- ap$corr
  AS$Ppaires[k] <- ap$po
  return(AS)
}

# AS$Crit <- rep(0, nc)      # max(abs(Corr))
# AS$Prob <- rep(0, nc)      # corrigé pour le nombre de corrélations
# AS$Ppaires <- rep(0, nc)   # les poids d'annulation
# AS$Corr <- matrix(0, AS$nv-2, nc)
# AS$doublet <- NULL


# CritLog <- function(px, G, combine, cible, K, to.opt = TRUE) {
#   pp <- c(exp(px), K)
#   SP <- G[, combine] %*% pp
#   SP <- SP %*% solve(sqrt(t(SP) %*% (SP)))
#   corr <- t(SP) %*% G
#   corr[setdiff(1:length(corr), setdiff(cible, combine))] <- 0
#   crit <- max(corr^2)
#   if (is.complex(crit)) {
#     crit <- 99e9
#   }
#   if(to.opt){
#     out <- crit
#   } else {
#     out <- list(crit = crit, corr = corr)
#   }
#   return(out)
# }
