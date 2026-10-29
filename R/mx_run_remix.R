mx_run_remix <- function(model, use_random_starts = TRUE) {
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
    if(use_random_starts){
      model_randomstarts <- try(random_starts(model, verbose = FALSE), silent = TRUE)

      if (mx_convergence(model_randomstarts, check = c("convergence", "fit"))$checks_pass) {
        model <- model_randomstarts
      } else {
        warning("Neither simulated annealing nor random starts converged.")
      }
    } else {
      warning("Simulated annealing failed. Consider using random starts.")
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
