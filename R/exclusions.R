exclusions <- function(AS,ote){
  if (!exists("exclusions",AS)) AS$exclusions <- list()
  out <- colonnesGS(AS,ote)
  # AS$excl_var <- excl <- out$ote
  # AS$excl_cibles <- cibl <- out$cibl
  AS$excl_var <- excl <- ote
  noms <- colnames(AS$GS)
  nc <- nchar(noms[length(noms)]) # nombre de lettres dans la dernière colonne
  cibl <- which(nchar(noms)==nc)       #(AS$nv+1):ncol(AS$GS)
  if (nc==1) cibl <- cibl[-ote]
  AS$excl_cibles <- cibl
  no <- length(ote)
  nv <- length(cibl)
  if (is.null(AS$satur_var)){
    satu <- matrix(NA,nrow=no,ncol=no)
    for (j in 1:(no-1))
      for (i in (j+1):no){
        sat <- asSatPaire(AS,ote[sort(c(i,j))])
        # sat <- tryCatch({
        # expr={asSatPaire(AS,ote[sort(c(i,j))])}
        # warning=function(w)
        #   browser()})
        satu[i,j] <- sat[1]
        satu[j,i] <- sat[2]
      }
    AS$satur_var <- sat <- rowSums(satu,na.rm=TRUE)/(no-1)
  }
  AA <- warnings();if(length(AA)>0) browser()
  for (j in 1:length(cibl)){
    covar <- covariables(AS,excl,cibl[j]) # seulement cov qui ne partagent pas de bruit
    if (!is.null(covar))
      AS <- exclusion(AS,c(cibl[j],covar))
  }
  BB <- warnings();if(length(BB)>0) browser()
  return(AS)
}