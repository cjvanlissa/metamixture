run_cached <- function(expr) {

  # Capture code without evaluating
  expr <- substitute(expr)

  # Get current chunk label
  label <- knitr::opts_current$get("label")

  if (is.null(label) || !nzchar(label)) {
    stop("run_cached() must be called inside a labeled R Markdown chunk.")
  }

  path <- file.path(cache_dir, paste0(label, ".rds"))

  if (isTRUE(run_everything)) {

    # Evaluate code in the calling environment
    result <- eval(expr, envir = parent.frame())

    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(result, path)

  } else {

    if (!file.exists(path)) {
      stop("Cached result not found: ", path)
    }

    result <- readRDS(path)
  }

  result
}
