mx_convergence <- function(model, check = c("convergence", "fit", "parameters", "standard_errors", "boundaries", "hessian"), ...) {

  if (!inherits(model, what = "MxModel")) {
    stop("'model' must be an OpenMx MxModel object.")
  }

  status <- model$output$status$code

  # OpenMx status codes 0 and 1 are generally acceptable.
  convergence_ok <- tryCatch(all(status %in% c(0L, 1L)), error = function(e)FALSE)

  # Fit value
  fit <- model$output$fit
  fit_ok <- length(fit) == 1L && is.finite(fit)

  # Parameter estimates
  pars <- OpenMx::omxGetParameters(model, free = TRUE)
  parameters_ok <- all(is.finite(pars))

  # Boundary estimates
  boundaries <- tryCatch(
    mx_boundary_estimates(model, ...),
    error = function(e) NULL
  )

  boundaries_ok <- !is.null(boundaries) &&
    !any(boundaries$boundary)

  # Standard errors, if available
  se <- tryCatch(
    summary(model)$parameters[, "Std.Error"],
    error = function(e) NULL
  )

  standard_errors_ok <- !is.null(se) &&
    length(se) > 0L &&
    all(is.finite(se))

  # Hessian diagnostics, if available
  hessian <- model$output$hessian

  hessian_ok <- NA

  if (!is.null(hessian) &&
      is.matrix(hessian) &&
      all(is.finite(hessian))) {

    eigenvalues <- eigen(
      (hessian + t(hessian)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    hessian_ok <- all(eigenvalues > 0)
  }

  # Overall convergence:
  checks_pass <- eval(parse(text = paste0(paste0(check, "_ok"), collapse = " && ")))

  # Additional diagnostics
  warnings <- character()

  if (!convergence_ok) {
    warnings <- c(warnings,
                  paste("Optimizer status:", status))
  }

  if (!fit_ok) {
    warnings <- c(warnings, "Non-finite fit value.")
  }

  if (!parameters_ok) {
    warnings <- c(warnings, "Non-finite parameter estimates.")
  }

  if (!boundaries_ok) {
    warnings <- c(warnings,
                  "One or more parameters are at their bounds.")
  }

  if (!standard_errors_ok) {
    warnings <- c(warnings,
                  "Standard errors are unavailable or non-finite.")
  }

  if (identical(hessian_ok, FALSE)) {
    warnings <- c(warnings,
                  "Hessian is not positive definite.")
  }

  list(
    checks_pass = checks_pass,
    convergence_ok = convergence_ok,
    fit_ok = fit_ok,
    parameters_ok = parameters_ok,
    boundaries_ok = boundaries_ok,
    standard_errors_ok = standard_errors_ok,
    hessian_ok = hessian_ok,
    warnings = warnings
  )
}
