mx_boundary_estimates <- function(fit, tol = 1e-6) {

  pars <- OpenMx::omxGetParameters(fit, free = TRUE)
  params <- fit@matrices

  out <- do.call(rbind, lapply(params, function(m) {
    if (!any(m@free)) return(NULL)
    idx <- which(m@free)
    data.frame(
      parameter = m@labels[idx],
      estimate = m@values[idx],
      lbound = m@lbound[idx],
      ubound = m@ubound[idx]
    )
  }))

  out$boundary <- (is.finite(out$lbound) &
    (out$estimate - out$lbound <= tol)) | (is.finite(out$ubound) &
    (out$ubound - out$estimate <= tol))
  rownames(out) <- NULL
  out
}
