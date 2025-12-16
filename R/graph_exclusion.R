graph_exclusion <- function(AAA,out=NULL){
  # AAA préparé par logit_params()
  # retourne un data.frame prêt pour imprimer une courbe logistique pour chaque condition
  if (is.null(out)){
    out <- data.frame(matrix(nrow=0,ncol=5))
    colnames(out) <- c("cond","B0","B1","rmin","rmax")
  }
  tmp <- data.frame(matrix(NA,nrow=1,ncol=5))
  colnames(tmp) <- c("cond","B0","B1","rmin","rmax")
#  nr <- nrow(out)
#  browser()
  for (k in 1:length(AAA)){
    tmp$cond <- names(AAA)[[k]]
    tmp$B0 <- AAA[[k]][1,1]
    tmp$B1 <- AAA[[k]][2,1]
    tmp$rmin <- AAA[[k]][1,5]
    tmp$rmax <- AAA[[k]][2,5]
#    nr <- nr+1
    out <- rbind(out,tmp)
#    curve(sigmoid(x,tmp$b0,tmp$B1),from=tmp$rmin, to=tmp$rmax )
  }
#  browser()
  # ggplot(out,x=c(rmin,rmax),aes(x=x))+
  #  stat_function(fun=sigmoid,args=list(B0,B1))
  return(out)
}