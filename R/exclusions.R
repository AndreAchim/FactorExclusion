exclusions <- function(AS,ote){
  if (!exists("exclusions",AS)) AS$exclusions <- list()
  out <- colonnesGS(AS,ote)
  AS$excl_var <- excl <- out$ote
  AS$excl_cibles <- cibl <- out$cibl
  no <- length(ote)
  if (ote[1] <= AS$nv) {
    satu <- matrix(NA,nrow=no,ncol=no)
    for (j in 1:(no-1))
      for (i in (j+1):no){
        sat <- asSatPaire(AS,ote[sort(c(i,j))])
        satu[i,j] <- sat[1]
        satu[j,i] <- sat[2]
      }
    AS$satur_var <- sat <- rowSums(satu,na.rm=TRUE)/(no-1)
  }
  con <- array(NA,dim=c(AS$nv,AS$nv-no,no))
  cor <- array(NA,dim=c(length(cibl),no-1,no))
  for (j in 1:length(cibl)){
    # nm1 <- noms[cibl[j]]
    covar <- covariables(AS,excl,cibl[j]) # seulement cov qui ne partagent pas de bruit
    # for (k in 1:no){
    #   nm <- paste0(nm1,noms[covar[k]])
  AS <- exclusion(AS,c(cibl[j],covar))
  }
  return(AS)
}