rapport_exclusion <- function(SF,noms){
  F <- SF$patron %*% SF$corrFct [,-1]
  q05 <- rep(NA,6)
  q01 <- rep(NA,6)
  #  for (Li in SF$paires){
  for (pa in 1:length(SF$paires)){
    Li <- SF$paires[[pa]]
    np <- ncol(Li)
    # attendu <- rep(NA,np)
    # for (k in 1:np){
    #   a <- Li[1,k]
    #   b <- Li[2,k]
    #   attendu[k] <- t(as.matrix(F[a,])) %*% as.matrix(F[b,])
    # }
    for (nom in noms) {
      print(nom)
      mt <- get(nom)[[pa]]
      for (k in 1:np){
        A <- mt[k,13:18,]
        for (j in 1:6){
          q05[j] <- mean(A[j,]<.05)
          q01[j] <- mean(A[j,]<.01)
        }
        if (q05[6]==.195) browser()
        B <- colMeans(A)
        m05 <- mean(B<.05)
        m01 <- mean(B<.01)
        # p05 <- rowMeans(mt[k,13:18,]<.05)
        # p01 <- rowMeans(mt[k,13:18,]<.01)
        C <- mt[k,1:6,]
        mCR <- rowMeans(C)
        C <- mt[k,7:12,]
        mr <- rowMeans(C)
        out <- sprintf('var: %d %d, moy(CR): %4.3f moy(r): %4.3f p<.05: %4.3f p<.01: %4.3f',Li[1,k],Li[2,k],mean(mCR),mean(mr),m05,m01)
        print(out)
        print(round(mCR,digits=3))
        print(round(mr,digits=3))
        print(q05)
        print(q01)
        # out <- sprintf('%d %d %4.3f  %4.3f %4.3f   %5.4f %5.4f',Li[1,k],Li[2,k],attendu[k],mCR,mr,p05,p01)
        # print(out)
      }
    }
  }
}
