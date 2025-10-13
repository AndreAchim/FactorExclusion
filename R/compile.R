compile <- function(NN){
  # NN est (k,22,brep) issu de proj_sim()
  d <- dim(NN)
  np <- d[1]  # nombre de paires
  nr <- d[3]
  inn <- 2.54/sqrt(nr)
  for (i in 1:np){
    co <- mean(NN[i,c(1,12),])
    m <- (NN[i,1,]+NN[i,12,])/2
    seo <- inn*sd(m)
    cv <- mean(NN[i,c(2,13),])
    m <- (NN[i,2,]+NN[i,13,])/2
    sev <- inn*sd(m)
    #    centres[,i] <- c(co,cv)
    centres <- c(co-seo,co,co+seo,cv-sev,cv,cv+sev)
    sig05 <- matrix(FALSE,nrow=nr,ncol=5)
    sig01 <- matrix(FALSE,nrow=nr,ncol=5)
    mo <- rep(NA,nr)
    for (k in 1:nr){
      p1 <- NN[i,7:11,k]
      p2 <- NN[i,18:22,k]
      pp <- (p1+p2)/2
      sig05[k,] <- pp<.05
      sig01[k,] <- pp<.01
    }
    print(round(c(i,centres),4))
    print(colSums(sig05))
    print(colSums(sig01))
    #  browser()    
    #       # Pour chaque test, récupère les rangs où p < .05 pour la condition k
    #       rangs_list <- list(rangs_test1[[k]], rangs_test2[[k]], rangs_test3[[k]], rangs_test4[[k]], rangs_test5[[k]])
    #       
    #       for (rangs in rangs_list) {
    #         if (length(rangs) > 0) {
    #           for (idx in rangs) {
    #             # idx est supposé être un vecteur d'indices pour les 5 premières dimensions
    #             arr[idx[1], idx[2], idx[3], idx[4], idx[5], k] <- arr[idx[1], idx[2], idx[3], idx[4], idx[5], k] + 1
    #           }
    #         }
    #       }
    #     }
    # }
    # centre[,i] <- c(mean(mo),median(mo))
  }
}