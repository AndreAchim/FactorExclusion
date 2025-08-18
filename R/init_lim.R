init_lim <- function(G,a,b){
  p <- .1*(-40.5:40.5)
  np <- length(p)
  pr_min <- matrix(0,np)
  pr_max <- matrix(0,np)
  for (k in 1:np){
      pro <- projections(p[k],G,a,b)
      pr_min[k] <- min(pro)
      pr_max[k] <- max(pro)
    }
#max_proj(p[k],G,a,b) * abs(p[k])  # pénaliser quand on s'éloigne de 0
  sm <- abs(pr_min + pr_max)
  sm[(pr_min * pr_max)>0] <- mean(sm)
  r_mi <- which.min(sm * abs(p))
  r_ <- which.min(sm)
 # browser()
#  if (r_mi<4)
    limi <- c(p[r_],p[r_mi])
#  else if (r_mi>np-4)
#    limi <- c(p[r_mi],200)
#    else
#      limi <- c(p[r_mi-3],p[r_mi+3])
  print(c(p[r_],r_,projections(p[r_],G,a,b)))
  return(list(limi=limi,proj_min_max=cbind(pr_min,pr_max)))
}