prep_FP <- function(SF){
  # SF est une liste avec <- les champs $patronm $corrFct et $paires. 
  # Ce dernier champ contient une liste de paires de rangs de variables (2,k)
  # tester en dehors du facteur exprimé dans les variables 1 et 2.
  # N est la taille des échantillons et nrep le nombre de jeux de données simulées
  # En sortie, on aura un dataframe
  N <- c(60,120,240,480,960,1920)
  nrep=1000
  FP <- matrix(0,nrow=3,ncol=6)
  if (isSymmetric(SF$corrFct))
    F <- SF$patron %*% sqrtm(SF$corrFct)
  else
    F <- SF$patron %*% SF$corrFct
  # réserver l'espace
  Li <- SF$paires[[1]]
  # effectuer les nrep simulations
  # options(warn=1)
  for (gr in 1:length(N)){
    set.seed(N[gr])
    id <- 100*N[gr] -1
    for (rp in 1:nrep){
      out <- gen_data(F,N[gr])
#      dt <- out$dt
      dt <- out  # parce que gen_data est devenu son ancienne version
      AS <- init_SCA(dt)
      AS <- exclusions(AS,c(1,2))
      for (p in 1:3){
        out <- test_proj(AS,Li[1,p],Li[2,p])
        p <- out$pr[2]
        q <- 0+(p<.05)
        if (!(q==0 | q==1)) browser()
        FP[p,rp]=FP[p,rp]+q
        li <- c(id+rp,gr,p,q)
        if (rp==1 && gr==1 && p==1){
          result <- data.frame(ID=li[1],Taille=li[2],Paire=li[3],FP=li[4])
        } else {
          result[nrow(result)+1,] <- li
        }
      }
    }
  }
  result$paire <- result$Paire-3.5
  result$Paire <- factor(result$Paire,ordered=FALSE)
  return(list(result,FP))
}
