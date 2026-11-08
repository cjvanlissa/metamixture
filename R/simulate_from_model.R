simulate_remix <- function(x, add_vi = TRUE){
  nsim = x$data$numObs
  classes <- num_classes(x)
  cp <- extract_postprob_remix(x)
  C <- apply(cp, 1, sample, x = seq_len(classes), size = 1L, replace = FALSE)
  do.call(rbind, lapply(unique(C), function(i){
    n_c <- sum(C == i)
    tau2 <- diag(x@matrices$S$values)
    if(add_vi){
      vis <- sample(x$V$values, size = n_c, replace = TRUE, prob = cp[, i])
    } else {
      vis = 0
    }
    data.frame(
      yi = stats::rnorm(n_c,
                 mean = x@matrices$M$values[1, i],
                 sd = sqrt(tau2[i] + vis)),
      vi = vis,
      C = i)
  }))
}

