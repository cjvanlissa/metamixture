mx_run_remix <- function(model) {
  # Simulated annealing
  model_anneal <- try(
    OpenMx::mxRun(
      OpenMx::mxModel(
        model,
        OpenMx::mxComputeSimAnnealing()
      ),
      intervals = FALSE,
      silent = TRUE
    ),
    silent = TRUE
  )

  # Finetune the results
  if (!inherits(model_anneal, "try-error")) {
    model_anneal@compute <- NULL
    model_anneal <- try(OpenMx::mxRun(model_anneal, intervals = FALSE, silent = TRUE), silent = TRUE)
  }

  # Random starts if necessary
  if (mx_convergence(model_anneal, check = c("convergence", "fit"))$checks_pass) {
    model <- model_anneal
  } else {
    model_randomstarts <- try(random_starts(model, verbose = FALSE), silent = TRUE)

    if (mx_convergence(model_randomstarts, check = c("convergence", "fit"))$checks_pass) {
      model <- model_randomstarts
    } else {
      stop("Neither simulated annealing nor random starts converged.")
    }
  }

  # Confidence intervals
  model@compute <- NULL
  result <- try(OpenMx::mxRun(model, intervals = TRUE, silent = TRUE), silent = TRUE)

  if (inherits(result, "try-error")) {
    warning("Confidence interval calculation failed.")
    return(model)
  }
  result
}
