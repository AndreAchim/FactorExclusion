rapport_exclusion_r <- function(SF,noms){
  F <- SF$patron %*% SF$corrFct [,-1]
  for (pa in 1:length(SF$paires)){
    Li <- SF$paires[[pa]]
    np <- ncol(Li)
    for (nom in noms) {
      print(nom)
      mt <- get(nom)[[pa]]
      N <- dim(mt)[3]
      df <- N-2
      for (k in 1:np){
        A <- colMeans(mt[k,7:12,])
        T <- A*sqrt(df/(1-A*A))
        P <- 2*pt(-abs(T),df)
        m05 <- N*mean(P<.05)
        m01 <- N*mean(P<.01)
        C <- mt[k,1:6,]
        mCR <- rowMeans(C)
        C <- mt[k,7:12,]
        mr <- rowMeans(C)
        out <- sprintf('var: %d %d, moy(CR): %4.3f moy(r): %4.3f Np<.05: %4.0f Np<.01: %4.0f',Li[1,k],Li[2,k],mean(mCR),mean(mr),m05,m01)
        print(out)
        # print(round(mCR,digits=3))
        # print(round(mr,digits=3))
        # print(q05)
        # print(q01)
        # out <- sprintf('%d %d %4.3f  %4.3f %4.3f   %5.4f %5.4f',Li[1,k],Li[2,k],attendu[k],mCR,mr,p05,p01)
        # print(out)
      }
    }
  }
}
