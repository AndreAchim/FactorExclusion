probNonDoublet <- function(AS, melange) {
  # retourne une probabilité nette (ajustée selon le nombre de corrélations
  # dont le maximum en valeurs absolues est utilisé ici)
  
  var <- setdiff(AS$pertinent, melange)
#  var <- AS$pertinent[-melange]  # dangereux si orpheline avant variable du mélange
  dl <- AS$N-2
  r <- max(abs(AS$R[var, melange])) # corrélation maximale des deux variables de mélange avec toutes les autres pertinentes
  t <- r/sqrt((1-r*r)/dl)
  pr <- 1 - (1-2*pt(-t,dl))^(2 * length(var))
  
  return(pr)
}
