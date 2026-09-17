seq_k_dim <- function(AS,tuple){
  ASk <- AS  # est-ce requis pour ne pas modifier AS en dehors de cette fonction?
  ASk$GS[,AS$orphelines] <- AS$colOrphelines # remettre dans $GS les colonnes des orphelines
  ASk$pertinent <- 1:ASk$nv
  if (length(tuple)==1)
    return(test_k1(AS))
  k <- length(tuple)
  fa <- fareg(ASk$R,k)
  Rreduit <- ASk$R
  diag(Rreduit) <- fa$h2
  var <- setdiff(1:ASk$nv,tuple)
  satur <- matrix(0,nrow=ASk$nv,ncol=k) # faut-il un ajustement pour les orphelines?
  dSat <- sqrt(fa$h2[tuple])
  diag(satur[tuple,]) <- dSat
  out <- NULL
  for (k in var){
    variables <- c(k,tuple)
    ou <- optim_tuple(ASk,variables)
    satur[k,] <- -ou$meilPoids[-1] * dSat
    out <- rbind(out,c(variables,ou$prob))
  }
  R <- Rreduit[tuple,tuple]
  d <- diag(1/sqrt(diag(R)))
  R <- d %*% R %*% d
  satur <- ajustePolarites(satur,R,Rreduit)
  pr <- min(out[,ncol(out)] * (nv-length(tuple)))
  # z <- qnorm(out[,ncol(out)])
  # pr <- pnorm(sum(z)/sqrt(length(z)));
  prochain <- out[which.min(out[,ncol(out)]),1]
  return(list(prob=pr,out=out,meil=tuple,prochain=prochain,satur=satur,Rfct=R))
}

ajustePolarites <- function(satur,Rfct,Rreduit){
  v <- which(rowSums(satur==0)==0)
  s <- polarites_de_correlations(Rreduit[v,v])
  for (k in 1:(10*length(v))){
    V <- satur %*% Rfct %*% t(satur)
    sv <- polarites_de_correlations(V[v,v])
    i <- which(sv != s)
    # sR <- sign(V[v,v])
    # aa <- rowSums(sR * abs(V[v])^2)
    # browser()
    if (length(i)==0)
      break
    else
      satur[v[i[1]],] <- -satur[v[i[1]],]
  }
  return(satur)  
}

polarites_de_correlations <- function(R,e=2){
  diag(R) <- 0
  sR <- sign(R)
  s <- sign(rowSums(sR * abs(R)^e))
}
