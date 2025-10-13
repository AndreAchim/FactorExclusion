exclusions <- function(AS,ote){
  AS$excl_var <- as.vector(ote)
  no <- length(ote)
  cibl <- setdiff(AS$pertinent,ote)
  AS$excl_cibles <- cibl
  con <- array(NA,dim=c(AS$nv,AS$nv-no,no))
  cor <- array(NA,dim=c(AS$nv,no,no))
  n_o <- no:1
  for (j in length(cibl):1){
    out <- exclusion(AS,c(cibl[j],ote))
    for (k in n_o){
      AS$excl_weig[[k]][j] <- out[[k]]$po
      AS$excl_crit[[k]][j] <- out[[k]]$crit 
#      AS$excl_cont[[k]][j] <- out[[k]]$contrast 
#      AS$excl_corr[[k]][j] <- out[[k]]$corr
      con[,j,k] <- out[[k]]$contrast
      cor[,,k] <- out[[k]]$corr
    }
  }
  for (k in n_o){
    AS$excl_cont[[k]] <- con[,,k]
    AS$excl_corr[[k]] <- cor[,,k]
  }
#  browser()
  return(AS)
}