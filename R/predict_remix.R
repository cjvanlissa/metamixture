predict_class_remix <- function(object, newdata = NULL, ...){
  if(is.null(newdata)){
    return(max.col(extract_postprob(object), ties.method = "first"))
  }
  stopifnot("Argument `newdata` must be a `data.frame` with columns `yi` and `vi`." = (is.data.frame(newdata) && isTRUE(all(names(newdata) %in% c("yi", "vi")))))
  cp <- class_densities(model = object, yi = newdata$yi, vi = newdata$vi, normalize = TRUE)
  max.col(cp, ties.method = "first")
}
