exclusions <- function(AS,ote){
  AS$excl_var <- as.vector(ote)
  AS$excl_cibles <- cibl <- setdiff(AS$pertinent,ote)
  no <- length(ote)
  satu <- matrix(NA,nrow=no,ncol=no)
  for (j in 1:(no-1))
    for (i in (j+1):no){
      sat <- asSatPaire(AS,ote[sort(c(i,j))])
      satu[i,j] <- sat[1]
      satu[j,i] <- sat[2]
    }
  AS$satur_var <- sat <- rowSums(satu,na.rm=TRUE)/(no-1)
  con <- array(NA,dim=c(AS$nv,AS$nv-no,no))
  cor <- array(NA,dim=c(length(cibl),no-1,no))
  n_o <- no:1
  for (j in length(cibl):1){
    out <- exclusion(AS,c(cibl[j],ote))
    for (k in n_o){
      AS$excl_weig[[k]][j] <- out[[k]]$po
      AS$excl_crit[[k]][j] <- out[[k]]$crit 
      con[,j,k] <- out[[k]]$contrast
      cor[j,,k] <- out[[k]]$corr
    }
  }
#  browser()
  for (k in n_o){
    AS$excl_cont[[k]] <- con[,,k]
    AS$excl_corr[[k]] <- cor[,,k]
  }
  return(AS)
}