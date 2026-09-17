asNettoie <- function(Gr) {
  # Le dernier regroupement de Gr est considéré gardé
  # Tous les sous-ensembles qui précèdent sont enlevés
  m <- length(Gr)
  if (m>0){
    for (k in seq(m, 1)) {
      for (j in Gr[[k]]) {
        for (i in seq_len(k - 1)) {
          if (j %in% Gr[[i]]) {
            Gr[[i]] <- numeric()
          }
        }
      }
    }
  }
#  browser()
  Gr[sapply(Gr,length)==0] <- NULL
  reste <- unlist(Gr[sapply(Gr,length)==1] )
  Gr[sapply(Gr,length)==1] <- NULL
  # Réordonner les groupes selon le rang de leurs variables
  ng <- length(Gr)
  prem <- numeric(ng)
  for (k in seq_len(ng)) {
    Gr[[k]] <- sort(Gr[[k]])
    prem[k] <- Gr[[k]][1]
  }
  
  Gr <- Gr[order(prem)]
  return(list(Gr = Gr, reste = reste))
}
