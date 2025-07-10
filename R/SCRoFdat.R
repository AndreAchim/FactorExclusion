
SCRoFdat <- function(R, N = NULL, seuils = c(.001,.25)){
  # Preliminaries ####
  AS <- list(dat = R)
  if(is.null(N)) {AS$N = nrow(R)} else {AS$N = N}
  if(isSymmetric(as.matrix(R))) {AS$R <- R} else {AS$R <- cov(R)}
  AS$seuils <- sort(seuils)
  AS$et <- sqrt(diag(AS$R)) # AA: Ceci est destiné à pouvoir exprimer la solution factorielle en termes des variables d'origine.
  iet <-  1 / AS$et
  AS$R <- AS$R*(iet %*% t(iet))
  if(det(AS$R) < 0) stop("\nLa matrice de corrélation n'a pas un déterminant positif.\n")
  rNEST <- Rnest::nest(AS$R, n = AS$N)
  rPA <- Rnest::pa(rNEST)
  cat("NEST suggests", rNEST$nfactors, "factors.\n")
  cat("PA_median suggests", rPA$nfactors, "factors.\n")
  AS$minFct <- min(rNEST$nfactors[[1]],rPA$nfactors)
  AS$nv <- ncol(AS$R)
  AS$GS <- chol(AS$R)
  AS <- asOrphelines(AS)
  return (AS)
}
  
