minima_locaux <- function(vect,x=NA,plus=0){
# avec plus=0, les points aux deux bouts de vect ne seront pas comptés comme des minima possibles
# mais avec plus>0 ils pourront l'être
# Avec plus < 0 ils le seront sûrement (pas une bonne idée)
# x, si présent doit être de la même longueur que vect; c'est un marqueur associéà chaque rang de vect
  BB <- c(vect[1]+plus,vect,vect[length(vect)]+plus)
  CC <- sign(diff(BB))
  double <- which(CC==0)
  DD <- CC[-double]
  EE <- diff(DD)
  rg <- which(EE==2)
  if (length(rg)>0) {
    for (k in 1:length(rg))
      rg[k] <- rg[k]+sum(rg[k]>=(double-1))
    if (is.vector(x))
      x <- x[rg]
    if (min(vect[c(1,length(vect))]) < min(vect[rg])) { # le minimum trouvé n'est que local 
      # rg <- which(EE==3)
      # if (length(rg) == 0) rg <- 0
      return(list(rang=0))  # retourner $rang = 0 (E n'est jamais 3)
    }
    else
      return(list(rang=rg,x=x,val=vect[rg]))
  } else {
    return(list(rang=0))
  }
}