p_signif <- function(ST,pa,paires,noms){
#  out <- prep_signif(ST,pa,paires,noms[1])
#  dat <- out$dat
#  Nsig <- out$Nsig
  Nsig <- prep_signif(ST,pa,paires,noms[1])
  for (k in 2:length(noms)){
    out <- prep_signif(ST,pa,paires,noms[k])
#    dat <- rbind(dat,out$dat)
#    Nsig <- rbind(Nsig,out$Nsig)
    Nsig <- rbind(Nsig,out)
  }
  # dat$Paire <- factor(dat$Paire,ordered=FALSE)
  # dat$ID <- factor(dat$ID,ordered=FALSE)
#  return(list(dat=dat,Nsig=Nsig))
  colnames(Nsig) <- c("x0y1","x0y2","x1y0","x1y2","x2y0","x2y1")
  Nsig <- rbind(Nsig,colSums(Nsig))
  rownames(Nsig) <- c(noms,"Total")
  return(Nsig)
}