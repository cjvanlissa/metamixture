check_convergence_cpp <- function(x, check = c("convergence", "fit", "parameters", "standard_errors", "boundaries", "hessian"), tol = 1e-6, ...) {

  if (!inherits(x, "remix_cpp")) {
    stop("'x' must be a remix_cpp object.")
  }

  # ECM convergence
  ecm_ok <- isTRUE(x$ecm$converged)

  # Joint optimization
  joint_attempted <- isTRUE(x$joint$attempted)
  joint_accepted <- isTRUE(x$joint$accepted)

  joint_ok <- joint_attempted &&
    joint_accepted &&
    identical(as.integer(x$joint$convergence), 0L)

  # Overall optimizer convergence
  convergence_ok <- if (joint_accepted) {
    joint_ok
  } else {
    ecm_ok
  }

  # Log-likelihood
  fit_ok <- length(x$logLik) == 1L &&
    is.finite(x$logLik)

  # Parameter estimates
  pars <- c(x$mu, x$tau2, x$pi)
  pi_ok <- all(x$pi >= 0) &&
    abs(sum(x$pi) - 1) <= 1e-6
  parameters_ok <- all(is.finite(pars)) && pi_ok

  # Equality constraints are not treated as boundary problems.
  # Check heterogeneity and mixing proportions.
  tau2_lower <- 0
  tau2_upper <- diff(range(x$yi))^2

  boundary_tau2 <- x$tau2 <= tau2_lower + tol |
    x$tau2 >= tau2_upper - tol

  boundary_pi <- x$pi <= tol

  boundaries_ok <- !any(boundary_tau2 | boundary_pi)

  # Standard errors
  se <- c(x$se_mu, x$se_tau2, x$se_pi)

  standard_errors_ok <- length(se) == 3L * x$K &&
    all(is.finite(se))

  # Observed hessian (information is negative hessian)
  hessian_ok <- NA

  if (!is.null(x$information) &&
      is.matrix(x$information) &&
      all(is.finite(x$information))) {

    ev <- eigen(
      (x$information + t(x$information)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    hessian_ok <- all(ev > 0)
  }

  # Diagnostic messages
  warnings <- character()

  if (!ecm_ok) {
    warnings <- c(warnings,
                  "ECM did not meet its stopping criterion.")
  }

  if (joint_attempted && !joint_ok) {
    warnings <- c(warnings,
                  paste("Joint optimization:",
                        x$joint$message))
  }

  if (!fit_ok) {
    warnings <- c(warnings, "Non-finite log-likelihood.")
  }

  if (!parameters_ok) {
    warnings <- c(warnings,
                  "Non-finite parameter estimates.")
  }

  if (!pi_ok) {
    warnings <- c(warnings,
                  "Mixing proportions are invalid.")
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
                  "Observed hessian is not positive definite.")
  }

  # Include existing C++ diagnostics
  if (length(x$notes)) {
    warnings <- c(warnings, x$notes)
  }
  pass_check <- eval(parse(text = paste0(paste0(check, "_ok"), collapse = " && ")))

  list(
    checks_pass = pass_check,
    convergence_ok = convergence_ok,
    fit_ok = fit_ok,
    parameters_ok = parameters_ok,
    boundaries_ok = boundaries_ok,
    standard_errors_ok = standard_errors_ok,
    hessian_ok = hessian_ok,
    warnings = unique(warnings)
  )
}
