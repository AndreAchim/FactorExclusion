exclutSequentiels <- function(Z, ff) {
  # exclut de ff un regroupement qui présume d'un autre pas inclu
  # Z fut préparée par hclust()
  # ff contient les groupes de lignes dans Z correspondant à des sous-ensembles
  # d'agrégations considérées incertaines (.001 > p > .25)
  if(length(ff) > 2){
    d <- nrow(Z)
    sur <- 1:(d + ff[[1]])
    
    for (f in (length(ff) - 1):2) {  # premier et derniers sous-ensembles jamais en cause
      g <- ff[[f]]
      z <- g + d + 1
      for (k in seq_along(g)) {  # tester chaque élément de g (ou de z)
        gr <- c(sur, setdiff(z, z[k]))
        if (!any(Z[g[k], 1] %in% gr) || !any(Z[g[k], 2] %in% gr)) {
          ff[[f]] <- NULL
          break
        }
      }
    }
    
    ff <- ff[!sapply(ff, is.null)]
  }
  return(ff)
}
