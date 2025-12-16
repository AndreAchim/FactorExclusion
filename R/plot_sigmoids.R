plot_sigmoids <- function(AA){
  #curve(sigmoid(x,tmp$b0,tmp$B1),from=tmp$rmin, to=tmp$rmax )
  # for (k in 1:nrow(AA)){
  #   xtp <- data.frame(x=seq(AA[k,rmin],AA[k,rmax],length.out=100))
  #   
  # }
  # x <- seq(0,.4,length.out=100)
  # ggplot(A_A,aes(x=x,linetype=cond))+
  #  stat_function(fun=sigmoid,args=list(A_A$B0,A_A$B1))
  # pl <- ggplot()
  # for (k in 1:6) {  #nrow(AA)){
  #   pl <- pl + xlim(AA[k,4],AA[k,5]) +
  #     stat_function(fun=sigmoid,args=list(B0=AA[k,2],B1=AA[k,3]),linetype=k)
  # }
  # pl
  x <- seq(0,.4,length.out=200)
  for (k in 1:nrow(AA)){
    b <- which(x>AA[k,4])[1]-1
    h <- which(x<AA[k,5])
    h <- h[length(h)]+1
    pr <- sigmoid(x[b:h],AA[k,2],AA[k,3])
    if (k==1)
      plot(x[b:h],pr,type="l",xlim=(c(0,.4)),ylim=c(0,1))
    else
      lines(x[b:h],pr,type="l")
  }
}