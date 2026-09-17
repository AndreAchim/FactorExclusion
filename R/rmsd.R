rmsd <- function(Est,Pop,Symmetric = TRUE, IncludeDiag = FALSE){
  D <- Est-Pop
  if(IncludeDiag) diag(D) <- 0
  sqrt(sum(D^2)/length(D))
}