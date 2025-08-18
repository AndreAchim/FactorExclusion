minima_locaux <- function(vect,x=NA,plus=0){
# avec plus=0, les points aux deux bouts de vect ne seront pas comptés comme des minima possibles
# mais avec plus>0 ils pourront l'être
# Avec plus < 0 ils le seront sûrement (pas une bonne idée)
# x, si présent doit être de la même ongueur que vect; c'est un marqueur associéà chaque grang de vect
  B <- c(vect[1]+plus,vect,vect[length(vect)]+plus)
  C <- sign(diff(B))
  double <- which(C==0)
  D <- C[-double]
  E <- diff(D)
  rg <- which(E==2)
  if (length(rg)>0) {
    for (k in length(rg))
      rg[k] <- rg[k]+sum(rg[k]>=(double-1))
    if (is.vector(x))
      x <- x[rg]
    if (min(vect[c(1,length(vect))]) < min(vect[rg])) # le minimum trouvé n'est que local 
      return(list(rang=which(E==3)))  # retourner $rang de longueur 0 (E n'est jamais 3)
    else
      return(list(rang=rg,x=x,val=vect[rg]))
  } else {
    return(list(rang=rg))
  }
}