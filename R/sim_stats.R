sim_stats <- function(fe,seuil=.05,N=NULL){
  # fe est une structure produite par FE_methodes2
  # la structure de sortie a les champs $CR, $CP et $nOK
  if (is.null(N)) N <- fe$N
  # t_crit <- qt(1 - seuil / 2, df = N - 2)
  # r_crit <- t_crit / sqrt(t_crit^2 + N - 2)
  # OK <- (fe$dx[,1]>r_crit) & (fe$dx[,2]==0)
  OK <- which(fe$Heywood==FALSE)
  nOK <- length(OK)
  fe_stats <- function(dt,att,pr){
    att <- atanh(att)
    np <- length(att)
    st <- matrix(nrow=6,ncol=np)
    sg <- matrix(nrow=2,ncol=np)
    for (k in seq_len(np)){
      adt <- atanh(dt[,k])-att[k]
      tt <- t.test(adt)
      st[,k] <- c(tt$statistic,tt$parameter,tt$p.value,tanh(c(tt$estimate,tt$conf.int[1:2])))
      sg[1,k] <- mean(pr[,k]<.05, na.rm=TRUE)
      sg[2,k] <- mean(pr[,k]<.01, na.rm=TRUE)
    }
    return(list(st=st,sg=sg))
  }
  aCR <- fe_stats(fe$CR[OK,],fe$CRatt,fe$pCR[OK,])
  # aCP <- fe_stats(fe$CP[OK,],fe$CRatt[fe$CRatt==0],fe$pCR[OK,])  # pour en avoir le bon nombre
  # return(list(CR=aCR$st,sgCR=aCR$sg,CP=aCP$st,sgCP=aCP$sg,nok=nOK))
  return(list(CR=aCR$st,sgCR=aCR$sg,nOK=nOK))
}