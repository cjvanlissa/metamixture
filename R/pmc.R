pmc <- function(x, ..., reps = 100, ci = .95){
  if(inherits(x[[1]], "MxModel")) return(pmc_mx(x = x, reps = reps, ci = ci, ...))
}

#' @importFrom stats median
pmc_mx <- function(x, ..., reps = 100, ci = .95){
  if (length(x) < 2L) {
    stop("At least two models are required.")
  }
  #nams <- c("comparison", "null", "alt", "null_stat", "alt_stat", "lb_dif",
  #          "ub_dif", "sig")
  nams <- tryCatch(vapply(x, function(m) m@name, character(1)), error = function(e){paste0("model ", seq_len(length(x)))})
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


pmc_mx2 <- function(x, ..., reps = 100L, ci = 0.95,
                   nboot = 1000L) {

  if (length(x) < 2L) {
    stop("At least two models are required.")
  }

  if (length(reps) != 1L || !is.finite(reps) ||
      reps < 2L || reps != as.integer(reps)) {
    stop("'reps' must be an integer >= 2.")
  }

  if (length(ci) != 1L || !is.finite(ci) ||
      ci <= 0 || ci >= 1) {
    stop("'ci' must be between 0 and 1.")
  }

  if (length(nboot) != 1L || !is.finite(nboot) ||
      nboot < 2L || nboot != as.integer(nboot)) {
    stop("'nboot' must be an integer >= 2.")
  }

  K <- length(x)

  nams <- tryCatch(
    vapply(x, function(m) m@name, character(1)),
    error = function(e) paste0("model ", seq_len(K))
  )

  if (anyDuplicated(nams)) {
    nams <- make.unique(nams)
  }

  # Observed effect sizes
  yi <- as.numeric(x[[1]]$Y$values)

  # One simulated score
  one_rep <- function(model) {

    sim <- simulate_remix(
      model,
      add_vi = TRUE
    )

    # Density estimated from simulated observations,
    # evaluated at observed effects
    score <- surprisal(
      x = sim$yi,
      y = yi,
      ...
    )

    if (length(score) != 1L ||
        !is.finite(score)) {
      stop("Surprisal must return one finite value.")
    }

    as.numeric(score)
  }

  # Simulate scores for each model
  rep_stat <- vapply(seq_along(x), function(i) {

    model <- x[[i]]

    if (requireNamespace("future.apply", quietly = TRUE)) {

      scores <- future.apply::future_replicate(
        n = reps,
        expr = one_rep(model),
        future.seed = TRUE,
        simplify = FALSE
      )

      unlist(scores, use.names = FALSE)

    } else {

      replicate(reps, one_rep(model))

    }

  }, FUN.VALUE = numeric(reps))

  colnames(rep_stat) <- nams

  # Mean expected surprisal
  means <- colMeans(rep_stat)
  names(means) <- nams

  # Comparisons:
  # successive models and each model against the first
  comparisons <- rbind(
    data.frame(
      comparison = "dif_seq",
      null_idx = seq_len(K - 1L),
      alt_idx = seq.int(2L, K)
    ),
    data.frame(
      comparison = "dif_one",
      null_idx = rep(1L, K - 1L),
      alt_idx = seq.int(2L, K)
    )
  )

  # Bootstrap CI for difference in mean scores
  compare_scores <- function(a, b) {

    delta <- mean(b) - mean(a)

    boot <- replicate(nboot, {

      aa <- sample(a, length(a), replace = TRUE)
      bb <- sample(b, length(b), replace = TRUE)

      mean(bb) - mean(aa)
    })

    limits <- stats::quantile(
      boot,
      probs = c((1 - ci) / 2, 1 - (1 - ci) / 2),
      names = FALSE
    )

    c(
      difference = delta,
      lb = limits[1],
      ub = limits[2]
    )
  }

  comparisons_stats <- t(vapply(
    seq_len(nrow(comparisons)),
    function(j) {

      compare_scores(
        rep_stat[, comparisons$null_idx[j]],
        rep_stat[, comparisons$alt_idx[j]]
      )

    },
    FUN.VALUE = c(
      difference = 0,
      lb = 0,
      ub = 0
    )
  ))

  out <- data.frame(
    comparison = comparisons$comparison,
    null = nams[comparisons$null_idx],
    alt = nams[comparisons$alt_idx],
    null_stat = means[comparisons$null_idx],
    alt_stat = means[comparisons$alt_idx],
    difference = comparisons_stats[, "difference"],
    lb = comparisons_stats[, "lb"],
    ub = comparisons_stats[, "ub"],
    sig = ifelse(
      comparisons_stats[, "ub"] < 0,
      "*",
      ""
    ),
    row.names = NULL
  )

  class(out) <- c("pmc_df", class(out))
  attr(out, "stat") <- "mean simulated-data KDE surprisal"
  attr(out, "scores") <- rep_stat

  out
}
