comp_test_corr <- function(N,pcrit=.05){
  dl <- N-2
  tcrit <- qt(pcrit/2,dl)
  t2 <- tcrit*tcrit/dl
  rt <- sqrt(t2/(t2+1))
  zcrit <- abs(qnorm(pcrit/2))
  rz <- tanh(zcrit/sqrt(N-3))
  return(list(rt,rz))
}