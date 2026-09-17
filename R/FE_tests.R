FE_tests <- function(FE=NULL,initial=c(.4,.5),saut=c(.1,.1),N=-c(125,500,2000)){
  if (is.null(FE))
    FE <- matrix(c(.4,.5,.4,.74,.4 ,-.4,.5,.5,.8,.6,.46,0,.54,.73,
                   0, 0, 0,.27,.55,.5,0,0,.3,.35,.6,.3,.4,.3,
                   0, 0, 0,  0,0,0,.6,.4,.4,.3,.4,.4,-.3,.4),ncol=3) 
  if (N[1]<0){
    N <- abs(N)
    set.seed(N[1])
  }
  sat <- initial
  out <- list()
  wrn <- getOption("warn")
  options(warn=-1)
  library(MASS)
  library(expm)
  while (max(sat)<.95) {
    FE[1:2,1] <- sat
    for (ns in N){
      cat(c(sat,ns),"\n")
      fe     <- FE_methodes2(FE,ns,c(1,2))
      OK_idx <- which(!fe$Heywood)           # indices des simulations valides
      sigCR  <- fe$pCR[OK_idx, ] < .05      # logiques, nOK × 66
      out[[length(out)+1]] <- list(sat=sat, ns=ns, stt=sim_stats(fe),
                                   sigCR=sigCR, OK_idx=OK_idx)
    }
    sat <- sat+saut
  }
  options(warn=wrn)
  CRatt_all <- fe$CRatt
  zro <- which(abs(CRatt_all)<1e-6)
  out[[length(out)+1]] <- list(CRatt=CRatt_all,zro=zro)
  # Rendre efg disponible globalement avant de sourcer les figures
  assign("efg", out, envir = .GlobalEnv)
  # source("R/fig_biais.R")             # produit Figure1.pdf et Figure1.png
  # source("R/build_detection_data.R")  # construit df_det, m0q, m3q
  # source("R/fig_detection.R")         # produit Figure2.pdf et Figure2.png
  return(out)
}