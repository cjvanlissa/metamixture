if(FALSE){
   library(OpenMx)
library(tidySEM)
library(metamixture)

set.seed(1234)

# Simulate data
n <- 300
pi <- c(0.6, 0.4)

# True latent class membership
C <- sample(1:2, n, replace = TRUE, prob = pi)

# Meta-analytic effect sizes
vi <- runif(n, 0.01, 0.05)
mu <- c(-1, 1)
tau2 <- c(0.1, 0.1)

yi <- rnorm(n, mean = mu[C], sd = sqrt(vi + tau2[C]))

# Auxiliary outcome with class-specific means
aux_mu <- c(0, 2)
aux <- rnorm(n, mean = aux_mu[C], sd = 0.5)

# Step 1: Estimate random-effects mixture
fit <- remix(
  yi = yi,
  vi = vi,
  classes = 2,
  mu = "varying",
  tau2 = "equal"
)

# Step 2: Estimate distal outcome model
res <- two_step(
  x = fit,
  model = "aux ~ 1",
  data = data.frame(aux = aux)
)

# Inspect results
table_results(res)

# Compare estimated auxiliary means with true values
OpenMx::omxGetParameters(res, free = TRUE)
}

two_step_internal_remix <- function(x, model, data, ...) {
  K <- num_classes(x)
  n <- nrow(x$Y$values)

  if (nrow(data) != n) {
    stop("Auxiliary data must have the same number of observations.")
  }

  if (!all(paste0("LogD", seq_len(K)) %in% names(x@algebras))) {
    stop("Expected class-specific LogD algebras in remix model.")
  }

  # Add class-specific auxiliary RAM models

  for (k in seq_len(K)) {

    aux <- label_auxiliary_parameters(model, k)

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

  # Construct joint class-specific log densities
  terms <- vapply(seq_len(K), function(k) {
    paste0(
      "LogD", k,
      " + log(weights[1,", k, "])",
      " + log(aux", k, ".fitfunction)"
    )
  }, character(1))

  # Replace LogTerms
  x <- OpenMx::mxModel(
    x,
    OpenMx::mxAlgebraFromString(
      paste0(
        "cbind(",
        paste(terms, collapse = ", "),
        ")"
      ),
      name = "LogTerms"
    )
  )

  # log-sum-exp calculation
  previous <- terms[1]

  for (k in 2:K) {

    step_name <- paste0("MaxStep", k)

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

  x <- OpenMx::mxModel(
    x,

    OpenMx::mxAlgebraFromString(
      previous,
      name = "RowMax"
    ),

    OpenMx::mxAlgebra(
      RowMax %*% Ones,
      name = "RowMaxMatrix"
    ),

    OpenMx::mxAlgebra(
      exp(LogTerms - RowMaxMatrix),
      name = "ExpTerms"
    ),

    OpenMx::mxAlgebra(
      ExpTerms %*% OnesCol,
      name = "RowSumExp"
    ),

    OpenMx::mxAlgebra(
      -2 * sum(RowMax + log(RowSumExp)),
      name = "Minus2LL"
    )
  )

  # Estimate auxiliary parameters
  x <- OpenMx::mxModel(
    x,
    OpenMx::mxFitFunctionAlgebra("Minus2LL")
  )

  OpenMx::mxRun(
    x,
    intervals = FALSE,
    silent = TRUE,
    ...
  )
}
