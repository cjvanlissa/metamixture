# set.seed(123)
# n = 1000
# k = 1000
# mu <- c(-.5, .5)
# tau2 <- c(.2, .8)
# c <- sample.int(2, size = k, replace = TRUE)
# Ti <- rnorm(k, mu[c], sd = sqrt(tau2[c]))
# df <- data.frame(t(sapply(Ti, function(ti){
#   x <- rnorm(n, ti)
#   c(yi = mean(x), vi = sd(x)/sqrt(n))
# })))
# library(OpenMx)
# library(tidySEM)
# metamixture:::remix(df$yi, df$vi, 1) -> res1
# metamixture:::remix(df$yi, df$vi, 2) -> res2
# metamixture:::remix(df$yi, df$vi, 3) -> res3
# res <- metamixture:::remix(df$yi, df$vi, 1:2)
# tmp <- summary(res[[2]])
# tmp$CI
# saveRDS(res, "res.rds")
# table_fit(res1)
# table_fit(res2)
# table_fit(res3)
#
# random_starts(fit2)

random_starts <- function(x,
                          nstarts = 20,
                          scale_A = c(-.5, .5),
                          scale_mat_dev = c(.001, 1.5),
                          jitter_thresholds = .25,
                          tolerance_minus2ll = 1e-6,
                          verbose = TRUE,
                          run = TRUE
) {

  if(is.null(.Random.seed)){
    warning("First initialize the random number generator")
    invisible(stats::rnorm(1))
  }
  out <- seed_values <- vector("list", length = nstarts)
  for(i in seq_len(nstarts)){
    seed_values[[i]] <- .Random.seed
    out[[i]] <- random_starts_remix_internal(x)
  }
  if(!run) return(out)
  out <- lapply(out, OpenMx::mxRun)
  tab_ll <- data.frame(seed = seq_along(seed_values), Minus2LogLikelihood = sapply(out, function(x){x$output$Minus2LogLikelihood}))
  if(sum(tab_ll$Minus2LogLikelihood - min(tab_ll$Minus2LogLikelihood) < tolerance_minus2ll) < 2){
    warning("The smallest Minus2LogLikelihood was not replicated.")
  }
  if(verbose) print(tab_ll)
  out <- out[[which.min(tab_ll$Minus2LogLikelihood)]]
  attr(out, "tab_ll") <- tab_ll
  attr(out, "seeds") <- seed_values
  return(out)
}

random_starts_remix_internal <- function(x) {
  # Set default values
  yi <- x$Y$values
  vi <- x$V$values
  n <- length(yi)
  dl_est <- dl_estimator(yi, vi)
  xtau2 <- dl_est["tau2"]
  # Plug in real values for numeric variables
  xmin <- min(x$Y$values)
  xmax <- max(x$Y$values)
  xrng <- xmax - xmin
  xrng[!is.finite(xrng) | xrng == 0] <- 1
  if(!is.finite(xtau2) | xtau2 == 0) xtau2 <- (xrng / 4)^2
  classes <- num_classes(x)
  for(i in seq_len(classes)){
    x@matrices[[paste0("Tau", i)]]$values[1,1] <- runif(1, min = 0.5 * xtau2, max = 1.0 * xtau2)
    x@matrices[[paste0("Mu", i)]]$values[1,1] <- runif(n = 1, min = xmin, max = xmax)
  }
  x@matrices$Eta$values[1, ] <- p_to_eta(p_from_parameters(x))
  x <- OpenMx::omxAssignFirstParameters(x)
  return(x)
}

# Calculate probabilities from assigned model parameters
p_from_parameters <- function(x){
  yi <- x$Y$values
  vi <- x$V$values
  classes <- num_classes(x)
  mu <- vapply(seq_len(classes), function(k) {
    as.numeric(x@matrices[[paste0("Mu", k)]]@values)
  }, numeric(1))

  tau2 <- vapply(seq_len(classes), function(k) {
    as.numeric(x@matrices[[paste0("Tau", k)]]@values)
  }, numeric(1))

  # Log of weighted component densities
  logdens <- sapply(seq_len(classes), function(k){
    dnorm(
      yi,
      mean = mu[k],
      sd = sqrt(vi + tau2[k]),
      log = TRUE
    )
  })
  row_max <- apply(logdens, 1, max)
  posteriors <- exp(logdens - row_max)
  return(colMeans(posteriors/rowSums(posteriors)))
}
