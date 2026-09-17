test_k1 <- function(AS){
  v <- which.min(AS$Prob)
  cp <- AS$Cpaires[,v]
  # a <- asPaireAvec(AS,cp[1])
  # b <- asPaireAvec(AS,cp[2])
  # pr <- c(sum(qnorm(AS$Prob[a])),sum(qnorm(AS$Prob[b])))
  ab <- rbind(asPaireAvec(AS,cp[1]),asPaireAvec(AS,cp[2]))
  pr <- rowSums(matrix(qnorm(AS$Prob[ab]),nrow=2))
  v <- which.min(pr)   # la variable directon retenue
  pr <- pnorm(pr[v])/sqrt(AS$nv-1)
  out <- matrix(c(ab[v,],AS$Prob[ab[v,]]),ncol=2)
  satur <- sat_sur_facteur(AS$R,AS$N,cp)$sat
  pa <- out[which.min(out[,2]),1]
  prochain <- setdiff(AS$Cpaires[,pa],cp[v])
  return(list(prob=pr,
              out=out,
              meil=cp[v],
              prochain=prochain,
              satur=satur,
              Rfct = 1)
           )
}


