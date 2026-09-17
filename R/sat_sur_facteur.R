sat_sur_facteur <- function(R, N, AB, seuil = 0.05) {
  # Estime les saturations factorielles de toutes les variables à partir
  # d'une paire (A, B) supposée exclusive au même facteur.
  #
  # Arguments :
  #   R     : matrice de corrélations (nv × nv)
  #   N     : taille d'échantillon
  #   AB    : vecteur de 2 indices : A = AB[1], B = AB[2]
  #   seuil : seuil de significativité pour le filtre des dénominateurs
  #
  # Retourne : list(sat, a_hat, b_hat, a_vec, b_vec, r_crit, n_valid_a, n_valid_b, details)
  #   sat    : vecteur de longueur nv — saturations estimées pour toutes les variables
  #            (NA pour A et B si on veut les distinguer — ici on y met a_hat et b_hat)
  #   details : data.frame des valeurs intermédiaires pour les variables hors AB
  ia     <- AB[1]
  ib     <- AB[2]
  ote <- AB
  nv     <- nrow(R)
  autres <- setdiff(seq_len(nv), AB)
  nva <- nv-2
  # R2p <- rep(NA,nva)
  # R2r <- rep(NA,nva)
  # R2s <- rep(NA,nva)
  # for (k in autres){
  #   rrr <- ratio_R2_stepwise(R,k,AB,N)
  #   R2p[k-2] <- rrr$R2_pair
  #   R2r[k-2] <- rrr$ratio
  #   R2s[k-2] <- rrr$R2_stepwise
  # }
  # Ri <- solve(R)
  # R2m <- 1 - 1/diag(Ri)
  # R2m[ote]
  # R2m <- R2m[-ote]

  # Valeur critique de |r| — calculée une seule fois pour ce N
  t_crit <- qt(1 - seuil / 2, df = N - 2)
  r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  r_AB <- R[ia, ib]
  r_AX <- R[ia, autres]
  r_BX <- R[ib, autres]
  
  # pour chaque corrélation de a et de B avec les autres variables
  sig_AX <- abs(r_AX) >= r_crit
  sig_BX <- abs(r_BX) >= r_crit
  
  # Condition de signe : r(A,B)*r(A,X)*r(B,X) > 0 garantit a²>0 et b²>0
  signe_ok <- r_AB * r_AX * r_BX > 0
  
  # x² = r(A,X)*r(B,X) / r(A,B)
  # signe de x : sign(r(A,X)) * sign(a_hat)
  # x estimé à 0 si :
  #   - aucune corrélation significative, OU
  #   - une seule significative mais r(A,X) et r(B,X) de signes opposés
  #     (contredit l'hypothèse d'un signal de X sur le facteur)
  #   - a_hat * b_hat == 0 (pas de facteur commun extractible — évite division par zéro)
  zero_x <- r_AX * r_BX * r_AB < 0
  # zero_x <- (!sig_AX & !sig_BX) | (xor(sig_AX, sig_BX) & r_AX * r_BX * r_AB < 0)
  x2_est <- r_AX * r_BX / r_AB
  x_est  <- ifelse(zero_x, 0.001,
                   sqrt(pmax(x2_est, 0)) * sign(r_AX) * sign(r_AB))
  # Vecteur complet des saturations (toutes variables)
  # sat        <- numeric(nv)
  # sat[ia]    <- a_hat
  # sat[ib]    <- b_hat
  # sat[autres] <- x_est
  # sat <- rbind(x_est,R2p,R2r)

  details <- data.frame(
    r_AX   = round(r_AX,  4),
    r_BX   = round(r_BX,  4),
    sig_AX = sig_AX,
    sig_BX = sig_BX,
    x_zero = zero_x,
    x_est  = round(x_est, 4)
  )
  # sat2 <- sat[1,]^2
  # R2p <- 1-(1-R2p)/2
  # R2r <- 1-(1-R2r)/2
  # R2s <- 1-(1-R2s)/2
  # Heywood <- c(any(sat2 > 1), any(sat2 > R2p),any(sat2 > R2r),any(sat2 > R2s),any(sat2 > R2m))
# if (Heywood[5]) browser()
  
  list(
    # sat       = sat,
    # r_crit    = r_crit,
    details   = details,
    Heywood   = any(abs(x_est) > 1)
  )
}
