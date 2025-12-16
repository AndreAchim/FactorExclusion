sigmoid <- function(r,B0,B1){
  1/(1+exp(-B1*r-B0))
}