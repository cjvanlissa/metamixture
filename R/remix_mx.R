remix_mx <- function(yi, vi,
                     classes = 2L,
                     mu = c("varying", "equal"),
                     tau2 = c("varying", "equal"),
                     run = TRUE) {
  intervals = TRUE

  if (!requireNamespace("OpenMx", quietly = TRUE)) {
    stop("Package 'OpenMx' is required.")
  }
  s <- search()
  if(!any(grepl("OpenMx", s, fixed = TRUE))){
    message("Please run `library(OpenMx) before running `remix()`.")
    return(NULL)
  }
  if(length(classes) > 1){
    out <- lapply(classes, function(i){
      remix_mx(yi = yi, vi = vi, classes = i, run = run)
    })
    attr(out, "tidySEM") <- "list"
    class(out) <- c("mixture_list", class(out))
    return(out)
  } else {
    # Input checks
    yi <- as.numeric(yi)
    vi <- as.numeric(vi)

    if (length(yi) != length(vi) ||
        length(yi) < 2L ||
        any(!is.finite(yi)) ||
        any(!is.finite(vi)) ||
        any(vi <= 0)) {
      stop("yi and vi must have equal lengths, non-missing, finite values, and vi > 0.")
    }

    if (length(classes) != 1L ||
        !is.finite(classes) ||
        classes < 1 ||
        classes != as.integer(classes)) {
      stop("'classes' must be a positive integer.")
    }

    classes <- as.integer(classes)
    n <- length(yi)

    # Starting values
    start_mu <- if (classes == 1L) {
      mean(yi)
    } else {
      as.numeric(stats::quantile(
        yi, probs = seq(0, 1, length.out = classes + 2L)[-c(1, classes + 2L)]
      ))
    }
    start_tau2 <- rep(max(0.001, stats::var(yi) - mean(vi)), classes)
    start_pi <- rep(1 / classes, classes)
    start_pi <- start_pi / sum(start_pi)

    # Observed data and known sampling variances
    objects <- list(
      OpenMx::mxMatrix(
        type = "Full", nrow = n, ncol = 1,
        free = FALSE, values = yi, name = "Y"
      ),
      OpenMx::mxMatrix(
        type = "Full", nrow = n, ncol = 1,
        free = FALSE, values = vi, name = "V"
      )
    )
    max_tau2 <- max(c(1e-7, diff(range(yi))^2))
    # Class-specific parameters and log-density algebras
    for (k in seq_len(classes)) {

      mu_name <- paste0("Mu", k)
      tau_name <- paste0("Tau", k)
      var_name <- paste0("S", k)
      logd_name <- paste0("LogD", k)

      objects <- c(objects, list(
        OpenMx::mxMatrix(
          type = "Full", nrow = 1, ncol = 1,
          free = TRUE, values = start_mu[k],
          lbound = min(yi),
          ubound = max(yi),
          labels = ifelse(mu == "varying", paste0("mu", k), "mu"),
          name = mu_name
        ),
        OpenMx::mxMatrix(
          type = "Full", nrow = 1, ncol = 1,
          free = TRUE, values = start_tau2[k],
          lbound = 1e-7,
          ubound = max_tau2,
          labels = ifelse(tau2 == "varying", paste0("tau2_", k), "tau2"),
          name = tau_name
        )
      ))

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          paste0("V + ", tau_name),
          name = var_name
        ),
        OpenMx::mxAlgebraFromString(
          paste0(
            "-0.5 * (log(2 * pi * ", var_name, ") + ",
            "(Y - ", mu_name, ")^2 / ", var_name, ")"
          ),
          name = logd_name
        )
      ))
    }

    # Mixing proportions and likelihood
    if (classes == 1L) {

      # No free mixing proportion
      objects <- c(objects, list(
        OpenMx::mxMatrix(
          type = "Full", nrow = 1, ncol = 1,
          free = FALSE, values = 1,
          name = "weights"
        ),
        OpenMx::mxAlgebraFromString(
          "-2 * sum(LogD1)",
          name = "Minus2LL"
        )
      ))

    } else {

      # classes-1 free logit parameters
      objects <- c(objects, list(
        OpenMx::mxMatrix(
          type = "Full",
          nrow = 1, ncol = classes - 1L,
          free = TRUE,
          values = log(start_pi[seq_len(classes - 1L)] / start_pi[classes]),
          labels = paste0("eta", seq_len(classes - 1L)),
          name = "Eta"
        )
      ))

      # Stable softmax: subtract maximum logit
      # The reference logit is zero.
      eta_expr <- paste0(
        "cbind(Eta, 0)"
      )

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          eta_expr, name = "Logits"
        ),
        OpenMx::mxAlgebraFromString(
          "exp(Logits - max(Logits)) / sum(exp(Logits - max(Logits)))",
          name = "weights"
        )
      ))

      # Log of weighted component densities
      # Build n x classes matrix of log(pi_k) + log density
      terms <- vapply(seq_len(classes), function(k) {
        paste0(
          "LogD", k, " + log(weights[1,", k, "])"
        )
      }, character(1))

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          paste0("cbind(", paste(terms, collapse = ", "), ")"),
          name = "LogTerms"
        )
      ))

      # Stable log-sum-exp over classes for each observation
      # Rowwise max via algebraic expression
      # Use exp(log density - rowwise maximum)
      # OpenMx's matrix algebra does not expose all R matrix
      # functions, so construct the expression explicitly.
      #
      # For each observation, calculate:
      # m_i = max_k LogTerms[i,k]
      # log f_i = m_i + log(sum_k exp(LogTerms[i,k]-m_i))
      #
      # An algebraic alternative is sequential logaddexp.
      lse <- terms[1]

      for (k in 2:classes) {
        a <- lse
        b <- terms[k]

        # log(exp(a) + exp(b)) =
        # max(a,b) + log(1 + exp(-abs(a-b)))
        #
        # Express max(a,b) elementwise using:
        # (a+b+abs(a-b))/2
        lse <- paste0(
          "((", a, " + ", b, " + abs(", a, " - ", b,
          "))/2 + log(1 + exp(-abs(", a, " - ", b, "))))"
        )
      }

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          paste0("-2 * sum(", lse, ")"),
          name = "Minus2LL"
        )
      ))
    }

    # Fit function
    objects <- c(objects, list(
      OpenMx::mxFitFunctionAlgebra("Minus2LL")
    ))

    # Confidence intervals
    if (intervals) {
      ci_names <- c(
        paste0("Mu", seq_len(classes)),
        paste0("Tau", seq_len(classes))
      )

      if (classes > 1L) {
        ci_names <- c(ci_names, "weights")
      }

      objects <- c(objects, list(
        OpenMx::mxCI(ci_names)
      ))
    }

    # Construct model
    model <- do.call(
      OpenMx::mxModel,
      c(list(paste0("RE_Mixture", classes)), objects, OpenMx::mxData(data.frame(yi), type = "raw"))
    )
    # Improve starting values
    if(classes > 1){
      strts <- try(remix_starts(model), silent = TRUE)
      if(!inherits(strts, what = "try-error")){
        model <- strts
      } else {
        message("Could not improve starting values.")
      }
    }

    if (run) {
      model <- OpenMx::mxModel(model, OpenMx::mxComputeSimAnnealing())
      res <- suppressWarnings(try(OpenMx::mxRun(model, silent = TRUE), silent = TRUE))
      if(!inherits(res, what = "try-error")){
        model <- res
        model@compute <- NULL
      }
      model <- OpenMx::mxRun(
        model, intervals = intervals
      )
    }
    attr(model, "tidySEM") <- "mixture"
    # OpenMx::mxSE("weights", fit_mx2)
    return(model)
  }
}
