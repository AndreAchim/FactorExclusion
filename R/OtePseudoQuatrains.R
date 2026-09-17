OtePseudoQuatrains <- function(AS) {
  brG <- 1
  while (brG <= length(AS$VG)) {
    copl <- AS$VG[[brG]]$coplan
    browser()
    if (nrow(copl) > 1) {
      coPlus <- classeTriplets(copl)
      codes <- unlist(unique(coPlus[, ncol(coPlus)]))
      rg <- which(codes %% 10 == 4)
      
      for (k in seq_along(rg)) { # Gestion possible de plusieurs plans de quatrains
        co <- coPlus[coPlus[, ncol(coPlus)] == codes[rg[k]], 2:(ncol(coPlus) - 1)]
        cens <- recens(co)
        
        if (any(cens[, 2] == 1)) {
          AS$VG[[brG]] <- NULL
          brG <- brG - 1
        }
      }
    }
    
    brG <- brG + 1
  }
  
  return(AS)
}