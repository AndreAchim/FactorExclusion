probNonDoublet <- function(AS, melange) {
  # Retourne une probabilité ajustée selon le nombre de corrélations.
  # Utilise le maximum en valeurs absolues des corrélations de la paire
  # avec toutes les autres variables pertinentes.
  var <- setdiff(AS$pertinent, melange)
  dl  <- AS$N - 2
  r   <- max(abs(AS$R[var, melange]))
  t   <- r / sqrt((1 - r * r) / dl)
  pr  <- 1 - (1 - 2 * pt(-t, dl))^(2 * length(var))
  return(pr)
}
