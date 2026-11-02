random_starts_internal_remix <- function(
    x,
    nstarts = 20,
    tolerance_minus2ll = 1e-6,
    verbose = TRUE,
    run = TRUE,
    ...
) {
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
  diag(x@matrices$S$values) <- stats::runif(n = length(diag(x@matrices$S$values)), min = 0.5 * xtau2, max = 1.0 * xtau2)
  x@matrices$M$values[1,] <- stats::runif(n = length(x@matrices$M$values[1,]), min = xmin, max = xmax)
  x@matrices$Eta$values[1, ] <- p_to_eta(p_from_parameters(x))
  x <- OpenMx::omxAssignFirstParameters(x)
  return(x)
}



# random_starts <- function(x,
#                           ...
# ) {
#
#   if(is.null(.Random.seed)){
#     warning("First initialize the random number generator")
#     invisible(stats::rnorm(1))
#   }
#   out <- seed_values <- vector("list", length = nstarts)
#   for(i in seq_len(nstarts)){
#     seed_values[[i]] <- .Random.seed
#     out[[i]] <- random_starts_internal(x, ...)
#   }
#   if(!run) return(out)
#   out <- lapply(out, OpenMx::mxRun, intervals = FALSE, silent = TRUE)
#   tab_ll <- data.frame(seed = seq_along(seed_values), Minus2LogLikelihood = sapply(out, function(x){x$output$Minus2LogLikelihood}))
#   if(sum(tab_ll$Minus2LogLikelihood - min(tab_ll$Minus2LogLikelihood) < tolerance_minus2ll) < 2){
#     warning("The smallest Minus2LogLikelihood was not replicated.")
#   }
#   if(verbose) print(tab_ll)
#   out <- out[[which.min(tab_ll$Minus2LogLikelihood)]]
#   attr(out, "tab_ll") <- tab_ll
#   attr(out, "seeds") <- seed_values
#   return(out)
# }
#
# random_starts_internal_remix <- function(x)
