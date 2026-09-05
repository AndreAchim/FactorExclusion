sat_sur_facteur <- function(R, N, AB, seuil = 0.05) {
  # Estime les saturations factorielles de toutes les variables à partir
  # d'une paire d'ancres (A, B) supposée exclusive au même facteur.
  #
  # Arguments :
  #   R     : matrice de corrélations (nv × nv)
  #   N     : taille d'échantillon
  #   AB    : vecteur de 2 indices — A = AB[1], B = AB[2]
  #   seuil : seuil de significativité pour le filtre des dénominateurs
  #
  # Retourne : list(details, Heywood)
  ia     <- AB[1]
  ib     <- AB[2]
  nv     <- nrow(R)
  autres <- setdiff(seq_len(nv), AB)

  t_crit <- qt(1 - seuil / 2, df = N - 2)
  r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  r_AB   <- R[ia, ib]
  r_AX   <- R[ia, autres]
  r_BX   <- R[ib, autres]

  sig_AX  <- abs(r_AX) >= r_crit
  sig_BX  <- abs(r_BX) >= r_crit
  signe_ok <- r_AB * r_AX * r_BX > 0
  zero_x   <- r_AX * r_BX * r_AB < 0
  x2_est   <- r_AX * r_BX / r_AB
  x_est    <- ifelse(zero_x, 0.001,
                     sqrt(pmax(x2_est, 0)) * sign(r_AX) * sign(r_AB))

  details <- data.frame(
    r_AX   = round(r_AX,  4),
    r_BX   = round(r_BX,  4),
    sig_AX = sig_AX,
    sig_BX = sig_BX,
    x_zero = zero_x,
    x_est  = round(x_est, 4)
  )

  list(
    details = details,
    Heywood = any(abs(x_est) > 1)
  )
}
