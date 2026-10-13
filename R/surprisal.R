surprisal <- function(x, y, aggregate = "sum", ...) {
  # Estimate density of x using KDE
  kde <- density(x, ...)

  # Evaluate estimated density at each y
  # rule = 2 extends boundary values outside KDE grid
  p <- approx(kde$x, kde$y, xout = y, rule = 2)$y

  # Calculate surprisal
  # Negative log kernel density
  out <- -log(pmax(p, .Machine$double.xmin))
  if(isFALSE(aggregate)){
    return(out)
  }
  if(identical(aggregate, "mean")) return(mean(out))
  return(sum(out))
}

pmc_remix <- function(x, ..., reps = 100, ci = .95){
  if (length(x) < 2L) {
    stop("At least two models are required.")
  }
  nams <- c("comparison", "null", "alt", "null_stat", "alt_stat", "lb_dif",
            "ub_dif", "sig")
  nams <- vapply(x_input, function(m) m@name, character(1))
  rep_stat <- vapply(seq_along(x_input), FUN.VALUE = numeric(reps), FUN = function(i) {
    if(requireNamespace("future.apply", quietly = TRUE)){
      do.call(c, future.apply::future_replicate(
        n = reps,
        future.seed = TRUE,
        simplify = FALSE,
        expr = {
          as.numeric(surprisal(
            x = drop(x_input[[1]]$Y$values),
            y = simulate_remix(x[[i]], add_vi = TRUE)$yi
          ))
        }
      ))
    } else {
      replicate(reps, {
        as.numeric(surprisal(
          x = drop(x_input[[1]]$Y$values),
          y = simulate_remix(x[[i]], add_vi = TRUE)$yi
        ))
      })
    }
  })
  # Get summary statistics and prepare output
  meds <- apply(rep_stat, 2, median)
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
  attr(out, "stat") <- which_stat
  return(out)
}
# library(OpenMx)
# library(tidySEM)
# df <- data_mix_ordinal
# df[1:4] <- lapply(df, ordered)
# mx_lca(data = df,
#        classes = 1:3) -> res
#
# pmc.mixture_list(res)
# data("empathy")
# df <- empathy[1:2]
# mx_profiles(data = df,
#             classes = 1:3) -> res
