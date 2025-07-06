
# SCRoF pairwise signal cancelleation
SCRoFpwsc <- function(AS){
  AS <- asPairesIndicatrices(AS)
  AS <- asGrappes(AS)
  if(is.null(AS$VG)) {cat('\nAucune annulation du signal par paire.\n')}
  AS <- asInitFct_Cor(AS) # CHECK 
}