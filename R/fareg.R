fareg <- function (R, numFactors = 1, facMethod = "rls") 
{
  if (!isSymmetric(R)) 
    stop("\nInput data must be a correlation matrix")
  if (!all.equal(as.vector(diag(R)), rep(1, dim(R)[1]))) {
    stop("\nInput data must be a correlation matrix")
  }
  PlausibleMethods <- c("rls", "rml")
  if (facMethod %in% PlausibleMethods == FALSE) {
    stop("\nInvalid argument for facMethod")
  }
  p <- nrow(R)
  communality <- "SMC"
  eigenVal <- eigen(R)$value
  if (min(eigenVal) <= 1e-06) {
    warning("Inverting the correlation matrix for SMC communality estimates requires a positive-definite matrix.")
    communality <- "maxr"
  }
  if (communality == "SMC") {
    Psi <- diag(diag(solve(R))^-1)
  }
  if (communality == "maxr") {
    R0 <- R
    diag(R0) <- 0
    h2temp <- apply(abs(R0), 2, max)
    Psi <- diag(1 - h2temp)
  }
  Lstart = 1
  FncRidgeULS <- function(L) {
    Ri <- R - (L * Psi)
    UDU <- eigen(Ri)
    U <- UDU$vectors
    D <- UDU$values
    D[D<0] <- 0   # ajout de AA pour éviter racine carrée d'une valeur négative
    if (numFactors == 1) {
      F <- U[, 1] * sqrt(D[1])
    }
    if (numFactors > 1) {
      F <- U[, 1:numFactors] %*% diag(sqrt(D[1:numFactors]))
    }
    Rhat <- F %*% t(F)
    rmsd(Ri, Rhat, Symmetric = TRUE, IncludeDiag = FALSE)
  }
  FncRidgeMLE <- function(L) {
    Ri <- R - (L * Psi)
    UDU <- eigen(Ri)
    U <- UDU$vectors
    D <- UDU$values
    if (numFactors == 1) {
      F <- U[, 1] * sqrt(D[1])
    }
    if (numFactors > 1) {
      F <- U[, 1:numFactors] %*% diag(sqrt(D[1:numFactors]))
    }
    Rhat <- F %*% t(F)
    diag(Rhat) <- 1
    RhatInv <- solve(Rhat)
    log(det(Rhat)) + sum(diag(R %*% RhatInv)) - log(det(R)) - 
      p
  }
  faSolution <- function(L) {
    Ri <- R - (L * Psi)
    UDU <- eigen(Ri)
    U <- UDU$vectors
    D <- UDU$values
    if (numFactors == 1) {
      F <- U[, 1] * sqrt(D[1])
    }
    if (numFactors > 1) {
      F <- U[, 1:numFactors] %*% diag(sqrt(D[1:numFactors]))
    }
    F
  }
  maxPsi <- max(diag(Psi))
  Lup <- 1/maxPsi + 0.001
  if (facMethod == "rls") {
    out <- optimize(f = FncRidgeULS, lower = 0.001, upper = Lup)
  }
  if (facMethod == "rml") {
    out <- optimize(f = FncRidgeMLE, lower = 0, upper = Lup)
  }
  Lopt <- out$minimum
  converged <- FALSE
  if (Lopt > 0 & Lopt < Lup) 
    converged <- TRUE
  loadings <- faSolution(Lopt)
  if (numFactors == 1) {
    loadings <- as.matrix(loadings)
    h2 <- as.vector(loadings^2)
  }
  if (numFactors > 1) {
    h2 <- apply(loadings^2, 1, sum)
  }
  Heywood <- FALSE
  if (max(h2) > 1) 
    Heywood <- TRUE
  list(loadings = loadings, h2 = h2, L = Lopt, converged = converged, 
       Heywood = Heywood)
}
