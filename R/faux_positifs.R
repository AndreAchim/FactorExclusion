faux_positifs <- function(DT,e=2,alpha=.05){
  fp <- c(0,0,0)
  mm <- c(1,1,1,1,1,1)
  rg <- list()
  for (k in 1:3){
    a <- colMeans(DT[[1]][k,3:4,] / DT[[1]][k,9:10,])^e  
    pr <- DT[[1]][k,13:18,]
    for (s in 1:length(a)){
      mm[3:4] <- a[s]
      if ((t(mm/sum(mm)) %*% pr[,s]) < alpha){
        browser()
        rg[length(rg)+1] <- s
        fp[k] <- fp[k] + 1
      }
    }
  }
  return(list(fp=fp,rg=rg))
}
