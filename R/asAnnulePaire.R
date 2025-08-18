asAnnulePaire <- function(AS, k) {
  # AS$GS(v,v) est la GrammSchmidt (en triangle supérieur) de R(v,v)
  # k est un scalaire
  # Optimise les variables de la rangée k de AS$Cpaires pour en minimiser le signal
  # Remplit AS$Crit[k], AS$Ppaires[k] et AS$Corr[,k]
  # à modifier: pour critère de proj maximale:
  # AS$Crit[k] est déjà le X2: N-1 fois la somme des corrélations au carré
  
  G <- AS$GS
  cible <- AS$pertinent
  melange <- AS$Cpaires[, k]
  lP <- 0
  pd <- probNonDoublet(AS, melange)
  
  if (pd < 0.05) {
    # r <- AS$R[melange[1], melange[2]]
    # K <- -sign(r)  # le poids fixe pour la 2e variable
    # lP <- 0  # logarithme de P=1
    # #POC: remove option and add them in optim
    # suppressWarnings({
    # lP <- optim(lP, fn = CritLog,
    #             method = "Nelder-Mead", 
    #             control =list(maxit = 1e9,
    #                           factr = 1e-6,
    #                           pgtol = 1e-6),
    #             G = G, combine = melange, K = K, cible = cible)$par
    # P <- exp(lP)
# AS$paires ne contient que des variables dans AS$pertinent
# Faut-il empêcher annule_paire() d'inclure les orpheliens parmi les témoins?
# Les 0 des orphelines ne vont jamais en faire la projection maximale
# mais on les aura en projections (nulles). Il ne faudait pas les compter
# dans le nombre de tests pour ajuter la probabilité
    ap <- annule_paire(G,melange[1],melange[2],AS$N)
    P <- ap$po
    } else {
    AS$doublet <- rbind(AS$doublet, melange)
    r <- AS$R[melange[1],melange[2]]
    t <- r/sqrt((1-r*r)/(AS$N-2))
    ap$prob <- 2*pt(-abs(t),AS$N-2)
    K <- 1
    }
if (length(AS$Corr[, k]) != length(ap$proj))  browser()
  
#  corr <- CritLog(lP, G, melange, cible, K, to.opt = FALSE)$corr
  AS$Prob[k] <- ap$prob
  AS$Crit[k] <- -log10(ap$prob)
  AS$Corr[, k] <- ap$proj
  AS$Ppaires[k] <- ap$po
  
  return(AS)
}

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
