simulate_remix <- function(x, add_vi = FALSE){
  nsim = x$data$numObs
  classes <- num_classes(x)
  cp <- tidySEM::class_prob(x, "individual")
  class_vi <- split(x$V$values, tidySEM::predict_class(x))
  C <- sample(seq_len(classes), size = nsim, replace = TRUE, prob = drop(x$weights$result))
  do.call(rbind, lapply(unique(C), function(i){
    n_c <- sum(C == i)
    tau2 <- x@matrices[[paste0("Tau", i)]]$values[1,1]
    if(add_vi){
      vis <- sample(x$V$values, size = n_c, replace = TRUE, prob = cp$individual[, i])
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
