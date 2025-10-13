prod_sim <- function(SF,N,nrep=100,ON=FALSE){
  # SF est une liste avec les champs $patronm $corrFct et $paires. 
  # Ce dernier champ contient une liste de paires de rangs de variables (2,k)
  # tester en dehors du facteur exprimé dans les variables 1 et 2.
  # N est la taille des échantillons et nrep le nombre de jeux de données simulées
  # En sortie, on aura une liste avec, pour chaque entrée de $paires, une matrice encore à déterminer
  if (N<0){
    N <- -N
    set.seed(N)
  }
  if (isSymmetric(SF$corrFct))
    F <- SF$patron %*% sqrtm(SF$corrFct)
  else
    F <- SF$patron %*% SF$corrFct
  # réserver l'espace
  result <- list()
  k <- 0
  for (Li in SF$paires){
    k <- k + 1
    np[k] <- ncol(Li)
    result[[k]] <- array(NA,c(np[k],22,nrep))
  }
  # effectuer les nrep simulations
  for (rp in 1:nrep){
#   print(rp)
    out <- gen_data(F,N,ON=ON)
    dt <- out$dt
    srce <- SC1(out$srce)
    AS <- init_SCA(dt)
    AS <- cree_zdat(AS)
    if (!is.null(SF$poids_attendus)) AS <- prepare_poids_fixes(AS,SF)
#    AS <- factor_exclusion(AS,c(1,2))
    AS <- exclusions(AS,c(1,2))
    k <- 0
#browser()
    for (Li in SF$paires){
      k <- k + 1
      for (p in 1:np[k]){
        out <- test_proj(AS,SF$paires[[k]][1,p],SF$paires[[k]][2,p],srce=srce)
        result[[k]][p,,rp] <- out
      }
    }
  }
  return(result)
}

prepare_poids_fixes <- function(AS,SF){
  v <- sort(unique(as.vector(SF$paires$zero))) +1  # +1 pour la colonne des rangs des covariables
  po <- SF$poids_attendus
  if ((length(v)+1) != (ncol(po))) browser()
  poids <- matrix(NA,nrow=nrow(po),ncol=AS$nv+1)
  poids[,1] <- po[,1]
  poids[,(v)] <- po[,-1]
  AS$poids_attendus <- poids
  return(AS)
}