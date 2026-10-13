simulate_remix <- function(x, add_vi = TRUE){
  UseMethod("simulate_remix", x)
}

#' @method simulate_remix MxModel
#' @export
simulate_remix.MxModel <- function(x, add_vi = TRUE){
  nsim = x$data$numObs
  classes <- num_classes(x)
  cp <- tidySEM::class_prob(x, "individual")[[1]][, seq_len(classes), drop = FALSE]
  C <- apply(cp, 1, sample, x = seq_len(classes), size = 1L, replace = FALSE)
  do.call(rbind, lapply(unique(C), function(i){
    n_c <- sum(C == i)
    tau2 <- x@matrices[[paste0("Tau", i)]]$values[1,1]
    if(add_vi){
      vis <- sample(x$V$values, size = n_c, replace = TRUE, prob = cp[, i])
      tau2 <- tau2 + vis
    } else {
      vis = 0
    }
    data.frame(
      yi = rnorm(n_c,
                 mean = x@matrices[[paste0("Mu", i)]]$values[1,1],
                 sd = sqrt(tau2)),
      vi = vis,
      C = i)
  }))
}

#' @method simulate_remix remix_cpp
#' @export
simulate_remix.remix_cpp <- function(x, add_vi = FALSE){
  nsim = x$n
  classes <- x$K
  cp <- x$posterior
  C <- apply(cp, 1, sample, x = seq_len(classes), size = 1L, replace = FALSE)
  do.call(rbind, lapply(unique(C), function(i){
    n_c <- sum(C == i)
    tau2 <- x$tau2[i]
    if(add_vi){
      vis <- sample(x$vi, size = n_c, replace = TRUE, prob = cp[, i])
      tau2 <- tau2 + vis
    } else {
      vis = 0
    }
    data.frame(
      yi = rnorm(n_c,
                 mean = x$mu[i],
                 sd = sqrt(tau2)),
      vi = vis,
      C = i)
  }))
}
