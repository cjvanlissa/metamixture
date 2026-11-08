auxiliary_internal_remix <- function(x, model, data, ...) {

  K <- num_classes(x)
  n <- nrow(data)

  # Check measurement data
  yi <- as.numeric(x$Y$values)
  vi <- as.numeric(x$V$values)

  if (length(yi) != n || length(vi) != n) {
    stop("Auxiliary data must match the measurement data.")
  }

  if (any(c("yi", "vi") %in% names(data))) {
    stop("Auxiliary variables cannot be named 'yi' or 'vi'.")
  }

  # Starting values based on the measurement-only model
  post_probs <- extract_postprob(x)

  # Add class-specific auxiliary models
  for (k in seq_len(K)) {

    aux <- label_auxiliary_parameters(model, k)

    aux <- mx_starts_auxiliary(
      model = aux,
      data = data,
      weights = post_probs[, k]
    )

    aux <- OpenMx::mxModel(
      aux,
      name = paste0("aux", k),
      OpenMx::mxData(
        observed = data,
        type = "raw"
      ),
      OpenMx::mxFitFunctionML(vector = TRUE)
    )

    x <- OpenMx::mxModel(x, aux)
  }

  # Combine measurement and auxiliary log densities.
  #
  # OpenMx vector ML fit functions return -2 log likelihood
  # contributions, so auxiliary log density = -0.5 * fit.
  terms <- vapply(seq_len(K), function(k) {
    paste0(
      "LogD", k,
      " + log(aux", k, ".fitfunction)",
      " + log(weights[1,", k, "])"
    )
  }, character(1))

  # Construct log-weighted component densities
  x <- OpenMx::mxModel(
    x,
    OpenMx::mxAlgebraFromString(
      paste0(
        "cbind(",
        paste(terms, collapse = ", "),
        ")"
      ),
      name = "AuxLogTerms"
    )
  )

  # Stable row-wise log-sum-exp
  previous <- terms[1]

  if (K > 1L) {
    for (k in seq.int(2L, K)) {

      step_name <- paste0("AuxMaxStep", k)

      expr <- paste0(
        "((", previous, ") + (", terms[k], ") + ",
        "abs((", previous, ") - (", terms[k], ")))/2"
      )

      x <- OpenMx::mxModel(
        x,
        OpenMx::mxAlgebraFromString(
          expr,
          name = step_name
        )
      )

      previous <- step_name
    }
  }

  x <- OpenMx::mxModel(
    x,

    OpenMx::mxAlgebraFromString(
      previous,
      name = "AuxRowMax"
    ),

    OpenMx::mxMatrix(
      type = "Unit",
      nrow = 1,
      ncol = K,
      name = "AuxOnes"
    ),

    OpenMx::mxAlgebra(
      AuxRowMax %*% AuxOnes,
      name = "AuxRowMaxMatrix"
    ),

    OpenMx::mxAlgebra(
      exp(AuxLogTerms - AuxRowMaxMatrix),
      name = "AuxExpTerms"
    ),

    OpenMx::mxMatrix(
      type = "Unit",
      nrow = K,
      ncol = 1,
      name = "AuxOnesCol"
    ),

    OpenMx::mxAlgebra(
      AuxExpTerms %*% AuxOnesCol,
      name = "AuxRowSumExp"
    ),

    OpenMx::mxAlgebra(
      -2 * sum(AuxRowMax + log(AuxRowSumExp)),
      name = "AuxMinus2LL"
    )
  )

  # Use the joint likelihood as the top-level objective.
  # Keep the original Minus2LL algebra unchanged.
  x <- OpenMx::mxModel(
    x,
    OpenMx::mxFitFunctionAlgebra("AuxMinus2LL")
  )

  # Reconcile shared auxiliary parameter labels
  x <- OpenMx::omxAssignFirstParameters(x)

  OpenMx::mxRun(
    x,
    intervals = FALSE,
    silent = TRUE,
    ...
  )
}
