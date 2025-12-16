linkage <- function(d, method = "complete", members = NULL){
#  rez <- hclust
  rez <- hclust(d, method, members)
  out <- rez$merge
#  out[out>0] <- out[out>0] + ncol(d)
  out[out>0] <- out[out>0] + ceiling(sqrt(2*length(d)))
  out <- abs(out)
  cbind(out, rez$height)
}
