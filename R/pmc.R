pmc <- function(x, ..., reps = 100, ci = .95){
  if(inherits(x[[1]], "MxModel")) return(pmc_mx(x = x, reps = reps, ci = ci, ...))
  pmc_cpp(x = x, reps = reps, ci = ci, ...)
}

#' @importFrom stats median
pmc_mx <- function(x, ..., reps = 100, ci = .95){
  if (length(x) < 2L) {
    stop("At least two models are required.")
  }
  #nams <- c("comparison", "null", "alt", "null_stat", "alt_stat", "lb_dif",
  #          "ub_dif", "sig")
  nams <- vapply(x, function(m) m@name, character(1))
  rep_stat <- vapply(seq_along(x), FUN.VALUE = numeric(reps), FUN = function(i) {
    if(requireNamespace("future.apply", quietly = TRUE)){
      do.call(c, future.apply::future_replicate(
        n = reps,
        future.seed = TRUE,
        simplify = FALSE,
        expr = {
          as.numeric(surprisal(
            x = drop(x[[1]]$Y$values),
            y = simulate_remix(x[[i]], add_vi = TRUE)$yi
          ))
        }
      ))
    } else {
      replicate(reps, {
        as.numeric(surprisal(
          x = drop(x[[1]]$Y$values),
          y = simulate_remix(x[[i]], add_vi = TRUE)$yi
        ))
      })
    }
  })
  # Get summary statistics and prepare output
  meds <- apply(rep_stat, 2, stats::median)
  names(meds) <- nams
  cis <- list(
    dif_seq = rep_stat[, 2:dim(rep_stat)[2], drop = FALSE] - rep_stat[, 1:(dim(rep_stat)-1L)[2], drop = FALSE],
    dif_one = rep_stat[, 2:dim(rep_stat)[2], drop = FALSE] - rep_stat[, rep(1, (dim(rep_stat)[2]-1L)), drop = FALSE]
  )
  cis <- do.call(rbind, lapply(cis, function(i){matrix(as.vector(apply(i, 2, stats::quantile, probs = c(((1-ci)/2), 1-((1-ci)/2)))), ncol = 2, byrow = TRUE)}))
  colnames(cis) <- c("lb", "ub")
  out <- data.frame(
    null = c(nams[1:(length(nams)-1L)], rep(nams[1], (length(nams)-1L))),
    alt = rep(nams[2:length(nams)], 2))
  out <- data.frame(comparison = rep(c("dif_seq", "dif_one"), each = nrow(out)/2),
                    out)
  out$null_stat <- meds[out$null]
  out$alt_stat <- meds[out$alt]
  out <- data.frame(out, cis, sig = c("", "*")[(apply(sign(cis), 1, sum) == -2)+1L])
  class(out) <- c("pmc_df", class(out))
  attr(out, "stat") <- "surprisal"
  return(out)
}

#' @importFrom stats median
pmc_cpp <- function(x, ..., reps = 100, ci = .95){
  if (length(x) < 2L) {
    stop("At least two models are required.")
  }
  #nams <- c("comparison", "null", "alt", "null_stat", "alt_stat", "lb_dif",
  #          "ub_dif", "sig")
  nams <- names(x)
  rep_stat <- vapply(seq_along(x), FUN.VALUE = numeric(reps), FUN = function(i) {
    if(requireNamespace("future.apply", quietly = TRUE)){
      do.call(c, future.apply::future_replicate(
        n = reps,
        future.seed = TRUE,
        simplify = FALSE,
        expr = {
          as.numeric(surprisal(
            x = drop(x[[1]]$yi),
            y = simulate_remix(x[[i]], add_vi = TRUE)$yi
          ))
        }
      ))
    } else {
      replicate(reps, {
        as.numeric(surprisal(
          x = drop(x[[1]]$yi),
          y = simulate_remix(x[[i]], add_vi = TRUE)$yi
        ))
      })
    }
  })
  # Get summary statistics and prepare output
  meds <- apply(rep_stat, 2, stats::median)
  names(meds) <- nams
  cis <- list(
    dif_seq = rep_stat[, 2:dim(rep_stat)[2], drop = FALSE] - rep_stat[, 1:(dim(rep_stat)-1L)[2], drop = FALSE],
    dif_one = rep_stat[, 2:dim(rep_stat)[2], drop = FALSE] - rep_stat[, rep(1, (dim(rep_stat)[2]-1L)), drop = FALSE]
  )
  cis <- do.call(rbind, lapply(cis, function(i){matrix(as.vector(apply(i, 2, stats::quantile, probs = c(((1-ci)/2), 1-((1-ci)/2)))), ncol = 2, byrow = TRUE)}))
  colnames(cis) <- c("lb", "ub")
  out <- data.frame(
    null = c(nams[1:(length(nams)-1L)], rep(nams[1], (length(nams)-1L))),
    alt = rep(nams[2:length(nams)], 2))
  out <- data.frame(comparison = rep(c("dif_seq", "dif_one"), each = nrow(out)/2),
                    out)
  out$null_stat <- meds[out$null]
  out$alt_stat <- meds[out$alt]
  out <- data.frame(out, cis, sig = c("", "*")[(apply(sign(cis), 1, sum) == -2)+1L])
  class(out) <- c("pmc_df", class(out))
  attr(out, "stat") <- "surprisal"
  return(out)
}
