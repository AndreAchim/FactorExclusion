asDistances <- function(AS){
  Dist <- matrix(0, AS$nv, AS$nv)
  for(k in 1:length(AS$Crit)){
    Dist[AS$Cpaires[2,k], AS$Cpaires[1,k]] <- AS$Crit[k]
  }
  # Dist[AS$orphelines,] <- NULL
  # Dist[,AS$orphelines] <- NULL
  if(length(AS$orphelines) > 0)  Dist <- Dist[-AS$orphelines, -AS$orphelines]
#  AS$Dist <- t(Dist * upper.tri(Dist))
#  AS$Dist <- Dist[Dist != 0]
#  AS$Dist <- as.dist(Dist+t(Dist), diag = FALSE, upper = FALSE)
  AS$Dist <- as.dist(Dist+t(Dist)) # retourne objet distances (triangle inférieur)
  # browser()
  # rez <- hclust(AS$Dist,"complete", NULL)
  # rez$merge
  return(AS)
}
