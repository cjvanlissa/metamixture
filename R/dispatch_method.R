dispatch_method <- function(x, generic, ...) {
  type <- attr(x, "tidySEM")
  fun <- switch(
    type,
    remix = paste0(generic, "_remix"),
    mixture = paste0(generic, "_mixture"),
    paste0(generic, "_MxModel")
  )

  get(fun, mode = "function")(x, ...)
}

table_fit <- function(x, ...) {

  if (methods::is(x, "MxModel")) {
    return(dispatch_method(x, "table_fit", ...))
  }
  # If it should handle other models too, then keep this to dispatch the regular
  # S3 methods:
  UseMethod("table_fit")
}
