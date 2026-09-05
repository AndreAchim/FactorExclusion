init_SCA <- function(R, N = NULL, seuils = NULL) {
  if (is.null(seuils)) seuils <- c(.001, .25)
  AS <- list(dat = R, seuils = seuils)
  if (is.null(N)) {AS$N <- nrow(R)} else {AS$N <- N}
  if (isSymmetric(as.matrix(R))) {AS$R <- R} else {AS$R <- cov(R)}
  AS$et  <- sqrt(diag(AS$R))
  iet    <- 1 / AS$et
  AS$R   <- AS$R * (iet %*% t(iet))
  if (det(AS$R) < 0)
    stop("\nLa matrice de corrélation n'a pas un déterminant positif.\n")
  if (is.null(colnames(AS$R)))
    rownames(AS$R) <- colnames(AS$R) <- letters[1:ncol(R)]
  AS$nv       <- ncol(AS$R)
  AS$pertinent <- 1:AS$nv
  AS$GS       <- chol(AS$R)
  colnames(AS$GS) <- letters[1:AS$nv]
  if (AS$nv > 26) colnames(AS$GS)[27:AS$nv] <- LETTERS[1:(AS$nv - 26)]
  AS$minFct   <- 1
  AS <- as_paires_indicatrices(AS)
  return(AS)
}
