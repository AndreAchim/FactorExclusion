test_prod_moy <- function(V1,V2=NULL){
  # Si V1 est présent seul, c'est un produit par cas
  # Si V1 et V2 présents, on fait leurs produits par cas
  if (!is.null(V2))
    V1 <- V1 * V2
  brut <- t.test(V1)
#  brut$kurtosis <- kurtosis(V1)-3
  tV1 <- sign(V1) * sqrt(abs(V1))
  trans <- t.test(tV1)
#  trans$kurtosis <- kurtosis(tV1)-3
#  return(list(brut=brut,trans=trans))
#  browser()
  N <- length(V1)
  be <- brut$estimate * N
  te <- trans$estimate * sqrt(N)
  bt <- brut$statistic
  tt <- trans$statistic
  bp <- brut$p.value
  tp <- trans$p.value
  p1b <- pt(brut$statistic,N-2)
  p2b <- pt(trans$statistic,N-2)
  out <- c(be,te,bt,tt,bp,tp,p1b,p2b)
  return(unname(out))
}