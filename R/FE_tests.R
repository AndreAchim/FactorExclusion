FE_tests <- function(FE = NULL, initial = c(.4, .5), saut = c(.1, .1),
                     N = -c(125, 500, 2000), seuil = .05) {
  # Simule les taux de faux positifs pour différentes combinaisons de
  # saturations d'ancres et de tailles d'échantillon.
  #
  # Arguments :
  #   FE      : matrice de saturations (défaut : matrice_EF intégrée au package)
  #   initial : saturations initiales des deux ancres
  #   saut    : incrément des saturations à chaque itération
  #   N       : tailles d'échantillon (valeurs négatives = graine fixe)
  #   seuil   : seuil de significativité pour la détection
  #
  # Retourne une liste dont chaque élément (sauf le dernier) contient :
  #   $sat    : saturations des ancres
  #   $ns     : N utilisé (négatif = graine fixe)
  #   $stt    : résultats de sim_stats()
  #   $sigCR  : matrice logique nOK × nc — détections à alpha = seuil
  #   $OK_idx : indices des simulations sans condition de Heywood
  # Le dernier élément contient $CRatt (population) et $zro (paires nulles).

  if (is.null(FE)) FE <- matrice_EF
  sat <- initial
  out <- list()
  wrn <- getOption("warn")
  options(warn = -1)

  while (max(sat) < .95) {
    FE[1:2, 1] <- sat
    for (ns in N) {
      cat(c(sat, ns), "\n")
      fe      <- FE_methodes2(FE, ns, c(1, 2))
      OK_idx  <- which(!fe$Heywood)
      sigCR   <- fe$pCR[OK_idx, ] < seuil
      out[[length(out) + 1]] <- list(
        sat    = sat,
        ns     = ns,
        stt    = sim_stats(fe),
        sigCR  = sigCR,
        OK_idx = OK_idx
      )
    }
    sat <- sat + saut
  }

  options(warn = wrn)
  CRatt_all <- fe$CRatt
  zro        <- which(abs(CRatt_all) < 1e-6)
  out[[length(out) + 1]] <- list(CRatt = CRatt_all, zro = zro)

  # Rendre efg disponible globalement pour les scripts de figures
  assign("efg", out, envir = .GlobalEnv)

  return(out)
}
