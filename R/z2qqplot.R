z2qqplot <- function(tronc,n,paire,E=1){
  nom <- paste(tronc,n,sep="")
  y <- get(nom)$z2[paire,]
  y <- qchisq(pchisq(y,1)^E,1)
  pr <- sprintf("%s, paire %d E=%1.2f",nom,paire,E)
## Q-Q plot for Chi^2 data against true theoretical distribution:
#qqplot(qchisq(ppoints(length(y)), df = 1), y, main = c(pr,expression("Q-Q plot for" ~~ {chi^2}[nu == 1])))
  qqplot(qchisq(ppoints(length(y)), df = 1), y, main = pr)
  lines(c(0,8),c(0,8),type="l")
# Add qq line
# library(ggplot2)
# ggplot2::last_plot() + qqline(y, distribution = function(p) qchisq(p, df = 1), prob = c(0.1, 0.6), col = 2)
#ggplot2::last_plot() + plot(c(0,8),c(0,8))
}
# my.qqline = qqline(y, distribution = function(p) qchisq(p, df = 3), prob = c(0.1, 0.6), col = 2)
# > my.qqline$data