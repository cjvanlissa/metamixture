remix_weights <- function(x){
  if(ncol(x$M$values) == 1L){
    return(1)
  }
  x$weights$result[1,]
}
