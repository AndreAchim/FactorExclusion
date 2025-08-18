max_abs <- function(x){
  if (length(x)<2)
    return(999)
  else
    max(abs(x))
}

sc1 <- function(x) {
# porte x à une somme de carrés de 1.0 sauf si tous des 0 dans x
  if (any(abs(x>1e-10)))
    x <- x / sqrt(sum(x^2))
  return(x)
  }

projections <- function(p,G,a,b){
  ctr <- sc1(p*G[a,]-G[b,]) %*% t(G[-c(a,b),])
  return(ctr)
}

max_proj <- function(p,G,a,b,fct=mamami){
  pro <- projections(p,G,a,b)
  return(fct(pro))
#  pro <- try(projections(p,G,a,b))
#  if(isTRUE(class(pro)=="try-error")) { browser() }
#  return(max(pro * pro))
#  return(mamami(pro))
#  return(abs(max(pro)+min(pro)))
}

paires_bifactorielles <- function(F,corr_crit=2){
# retourne toutes les paires de variables de F qui n'impliquent pas un seul facteur de F
# sauf celles dont le signal corrèle > corr_crit
  nv <- nrow(F)
  nf <- ncol(F)  # planterait si nf==1
  np <- nv
  for (k in 2:nf)
    np <- np + (nv-k)
#  paires <- matrix(0,np,2)
  prod_pa <- rep(0,np)
  cor_signal <- rep(0,np)
  p <- 0
  for (f1 in 1:(nf-1))
    for (j in 1:(nv-1))
      for (f2 in (f1+1):nf)
        for (k in (j+1):nv)
          if (F[j,f1]!=0 && F[k,f2]!=0) {
            pa <- nv*k+j
            if (p==0 || !any(prod_pa[1:p]==pa)){
              co <- sc1(F[j,]) %*% sc1(F[k,])
            if (any((F[j,] * F[k,]) != 0) && co < corr_crit) {
                p <- p+1
                prod_pa[p] <- pa
                cor_signal[p] <- co
              }
            }
          }
  paires <- matrix(0,p,3)
  paires[,1] <- prod_pa[1:p] %% nv
  paires[,2] <- (prod_pa[1:p]-paires[,1])/nv
  paires[,3] <- cor_signal[1:p]
  return(paires)
}
  
N_pour_contrastes <- function(F,R=NULL,p_crit=c(.25,.05,.001),corr_crit=.95){
# F(nv,nf) est matrice de patrons, R(nf,nf) les corrélatiosn factorielles
# lamatrice identité sera substituée si R n'est pas une matrice
# nf>1 est obligatoire
# utilisra G <- F %*% sqrtm(R)
# a et b les rangs des variables du contraste dans F (poids*G[a,]-G([b,]))
  if (!is.matrix(F)) stop("F doit être une matrice (nv,nf) avec nf > 1.")
  nv <- nrow(F)
  nf <- ncol(F)
  if (!is.matrix(R)) R <- diag(nf)
  G <- F %*% sqrtm(R)
  paires <- paires_bifactorielles(F,corr_crit)
  np <- nrow(paires)
  proj <- matrix(0,np,nv-2)
  crit <- rep(0,np)
  poids <- rep(0,np)
  for (k in 1:np){
    a <- paires[k,1]
    b <- paires[k,2]
    out <- optimize(max_proj,c(-20,20),G,a,b)
    crit[k] <-out$objective
    poids[k] <- out$minimum
    proj[k,] <- projections(poids[k],G,a,b)
  }
  pc <- (1-p_crit)^(1/(nv-2)) # 1-prob par test pour que tous les tests soient non significatifs
  X2_crit <- qchisq(pc,1)
  N <- ceiling(1+X2_crit / min(crit))
  return(list(N=N,paires=paires,crit=crit,poids=poids,proj=proj))
}