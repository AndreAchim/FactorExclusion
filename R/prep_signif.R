prep_signif <- function(ST,pa,rgs,nom,alpha=.05){
  # ST (peut être ST_Fl), pa (rang de ST$paires), rgs une liste de rangs dans cette matrice [k,20,nrep]
  paires <- ST$paires[pa]
  paires <- paires[[1]]     # [,rgs]
  mt <- get(nom[[1]])[pa][[1]]
  N <- dim(mt)[3]
  np <- length(rgs)
  Nsig <- rep(0,np)  # le décompte des positifs pour chaque paire dans la condition 'nom'
  id1 <- as.numeric(substring(nom,2))
  ta <- which(id1==c(60,120,240,480,960,1920))
  ID <- rep(100*id1+(0:(N-1)),each=np)
  Taille <- rep(ta-3.5,N*np)
  Paire <- rep(rgs,N)
  Sign <- rep(NA,N*np)
  dat <- data.frame(ID,Taille,Paire,Sign)
#  browser()
  dl <- id1-2
  sdl <- sqrt(id1-3)
  rng <- 0  # rangée dans dat
  for (k in 1:N){
    for (p in 1:np){
      rng=rng+1
      z=mean(atanh(mt[rgs[p],7:12,k]))*sdl
      pr <- 2*pnorm(-abs(z))
      si <- pr < alpha
      if (si) {
        Nsig[p] <- Nsig[p]+1
        }
      dat[rng,4] <- 0+si
    }
  }
  return(list(dat=dat,Nsig=Nsig))
}