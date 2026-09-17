optim_base <- function(AS, tuple, h2 = NULL) {
  # Teste si les k variables du tuple suffisent à expliquer toutes les autres.
  # Normalise les signaux à amplitude 1 (une seule fois), utilise dernier_poids
  # pour la contrainte de communauté, et agrège les probabilités par Stouffer.
  #
  # Retourne : list(prob, details)
  #   prob    : probabilité agrégée (Stouffer) pour H0 : k dimensions suffisent
  #   details : matrice indexée par variable (R2, prob, w1, ..., wk)
  #             poids rétablis à l'échelle originale; R2 et prob = NA pour les directions

  k <- length(tuple)
  methode = "Nelder-Mead"
  # Communalités et normalisation — calculés une seule fois
  if (is.null(h2)) h2 <- fareg(AS$R, k)$h2
  ampl <- sqrt(h2)
  GS_norm <- sweep(AS$GS, 2, ampl, "/")

  # Signaux-direction normalisés (base fixe pour toutes les variables cibles)
  DIRS <- GS_norm[, tuple, drop = FALSE]  # n × k

  var <- setdiff(1:AS$nv, tuple)   # rangs des variables à expliquer par DIRS
  nv  <- length(var)

  out <- matrix(NA, nrow = AS$nv, ncol = 2 + k,
                dimnames = list(NULL, c("R2", "prob",
                                        paste0("w", seq_len(k)))))
  out[tuple, 3:(2+k)] <- 0
  diag(out[tuple, 3:(2+k)]) <- ampl[tuple]

  for (v in var) {
    TARGET  <- GS_norm[, v]
    cible   <- AS$GS[, v]

    # Témoins : toutes les variables sauf la cible et les directions
    TEMOINS <- GS_norm[, setdiff(var, v), drop = FALSE]
    Rinv    <- solve(t(TEMOINS) %*% TEMOINS)

    crit <- function(p_free, DIRS, cible) {
      wk <- dernier_poids(p_free, DIRS, cible)
      if (is.na(wk)) return(1)
      poids    <- c(p_free, wk)
      contrast <- sc1(TARGET - DIRS %*% poids)
      z  <- t(contrast) %*% TEMOINS
      as.numeric(z %*% Rinv %*% t(z))
    }

    result <- optim(rep(0, k - 1), crit, gr = NULL, DIRS, cible,
                    method = methode,
                    control = list(parscale = rep(0.1, k - 1)))

    # Recalculer les poids finaux à partir de result$par
    wk         <- dernier_poids(result$par, DIRS, cible)
    poids_orig <- c(result$par, wk) * ampl[v]

    out[v, ] <- c(result$value,
                  prob_R2(result$value, ncol(TEMOINS), AS$N),
                  poids_orig)
  }

  # Agrégation de Stouffer
  z_comb <- sum(qnorm(out[var, "prob"])) / sqrt(nv)
  list(prob    = pnorm(z_comb),
       details = out)
}
