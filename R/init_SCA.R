init_SCA <- function(R, N = NULL){
  AS <- list(dat = R,seuils=c(.001,.25))
  if(is.null(N)) {AS$N = nrow(R)} else {AS$N = N}
  if(isSymmetric(as.matrix(R))) {AS$R <- R} else {AS$R <- cov(R)}
  AS$et <- sqrt(diag(AS$R)) # AA: Ceci est destiné à pouvoir exprimer la solution factorielle en termes des variables d'origine.
  iet <-  1 / AS$et
  AS$R <- AS$R*(iet %*% t(iet))
#  browser()
  if(det(AS$R) < 0) stop("\nLa matrice de corrélation n'a pas un déterminant positif.\n")
  AS$nv <- ncol(AS$R)
  AS$pertinent <- 1:AS$nv  # avant d'exclure les variables orphelines
  AS$GS <- chol(AS$R)
#  AS <- asOrphelines(AS) 
  # POC: retirer un conditionnel ici, utiliser deux fois ####
#  AS$pertinent <- setdiff(1:AS$nv, AS$orphelines)
  AS$minFct <- 1  # en attendant de restituer rNEST
# AS$N <- 2000 #### ATTENTION! à retirer ######
  AS <- asPairesIndicatrices(AS)
  return(AS)
}