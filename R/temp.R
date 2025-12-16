temp <- function(AS,cx,cy,x,y) {
AS <- cree_zdat(AS)
C1x <-  sc1(AS$zdat[,cx] - AS$excl_weig[[1]][x]*AS$zdat[,AS$excl_var[1]])
C1y <-  sc1(AS$zdat[,cy] - AS$excl_weig[[1]][y]*AS$zdat[,AS$excl_var[1]])
C2x <-  sc1(AS$zdat[,cx] - AS$excl_weig[[2]][x]*AS$zdat[,AS$excl_var[2]])
C2y <-  sc1(AS$zdat[,cy] - AS$excl_weig[[2]][y]*AS$zdat[,AS$excl_var[2]])
srce <-SC1(AS$srce)
browser()
cc <- c(1,2,4,5,10,11)
t(C1x) %*% srce[,cc]
t(C1y) %*% srce[,cc]
t(C2x) %*% srce[,cc]
t(C2y) %*% srce[,cc]
}