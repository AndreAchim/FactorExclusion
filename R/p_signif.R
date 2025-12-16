p_signif <- function(ST,pa,paires,noms){
  out <- prep_signif(ST,pa,paires,noms[1])
  dat <- out$dat
  Nsig <- out$Nsig
  for (k in 2:length(noms)){
    out <- prep_signif(ST,pa,paires,noms[k])
    dat <- rbind(dat,out$dat)
    Nsig <- rbind(Nsig,out$Nsig)
  }
  dat$Paire <- factor(dat$Paire,ordered=FALSE)
  dat$ID <- factor(dat$ID,ordered=FALSE)
  return(list(dat=dat,Nsig=Nsig))
}