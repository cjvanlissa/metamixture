remix_starts <- function(model,
                           splits = NULL,
                           ...){
  if(!isTRUE(requireNamespace("OpenMx", quietly = TRUE))) {
    return(NULL)
  }
  stopifnot("mxModel is not a `remix` model." = is_remix(model))

  mu_equal <- length(unique(model$M$labels[1,])) == 1L
  tau_equal <- length(unique(diag(model$S$labels))) == 1L
  dl_overall <- dl_estimator(y = model@matrices$Y$values, v = model@matrices$V$values)
  classes <- num_classes(model)
  if(classes == 1L){
    model@matrices$M$values[1,1] <- dl_overall["mu"]
    model@matrices$S$values[1,1] <- dl_overall["tau2"]
    return(model)
  }

  yi <- model@matrices$Y$values
  if(is.null(splits)){
    if (mu_equal && !tau_equal) {
      cluster_data <- (as.numeric(yi) - dl_overall["mu"])^2 / (as.numeric(model$V$values) + dl_overall["tau2"])
    } else {
      cluster_data <- yi
    }

    splits <- try({stats::kmeans(x = cluster_data, centers = classes)$cluster}, silent = TRUE)
    if(inherits(splits, "try-error")){
      message("Could not initialize clusters using K-means, switching to hierarchical clustering.")
      splits <- try({stats::cutree(stats::hclust(dist(yi)), k = classes)}, silent = TRUE)
      if(inherits(splits, "try-error")){
        stop("Could not initialize clusters using hierarchical clustering. Consider using a different clustering method, or imputing missing data.")
      }
    }
  }
  stopifnot("Number of unique values in splits must be identical to the number of latent classes." = length(unique(splits)) == classes)

  # Handle small cats
  splits <- repair_small_clusters(splits = splits, yi = cluster_data, min_size = 2L)

  if (classes != length(unique(splits))) {
    stop("Argument 'splits' does not identify a number of groups equal to 'classes'.")
  }
  # Order splits from largest number to smallest
  tab_split <- sort(table(splits), decreasing = TRUE)
  splits <- as.integer(ordered(splits, levels = names(tab_split)))

  dl_class <- lapply(seq_len(classes), function(i){
      dl_estimator(y = model@matrices$Y$values[splits == i], v = model@matrices$V$values[splits == i])
  })

  if(mu_equal){
    model$M$values[1,] <- dl_overall["mu"]

  } else {
    model$M$values[1,] <- sapply(dl_class, `[`, "mu")
  }

  if (tau_equal) {
    diag(model$S$values) <- rep(dl_overall["tau2"], classes)
  } else {
    diag(model$S$values) <- sapply(dl_class, `[`, "tau2")
  }

  p_start <- as.numeric(tab_split) / length(splits)
  model@matrices$Eta$values[1, ] <- p_to_eta(p_start)
  model <- OpenMx::omxAssignFirstParameters(model)
  return(model)
}



dl_estimator <- function (y, v)
{
  if(length(y) == 1L){
    return(c(mu = y, tau2 = 0))
  }
  Wi <- 1/v
  sum_wi <- sum(Wi)
  Q <- sum(Wi * y^2) - (sum(Wi * y)^2/sum_wi)
  df <- length(y) - 1
  C <- sum_wi - (sum(Wi^2)/sum_wi)
  tau2_est <- max(1e-3, (Q - df)/C)
  Wistar <- 1/(v + tau2_est)
  y_pooled <- sum(Wistar * y)/sum(Wistar)
  return(c(mu = y_pooled, tau2 = tau2_est))
}

p_to_eta <- function(p) {
  log(p[-length(p)] / p[length(p)])
}
eta_to_p <- function(eta) {
  logits <- c(eta, 0)
  p <- exp(logits - max(logits))
  p / sum(p)
}



repair_small_clusters <- function(splits, yi, min_size = 2L) {
  splits <- as.integer(factor(splits))
  K <- length(unique(splits))

  repeat {
    counts <- tabulate(splits, nbins = K)

    small <- which(counts < min_size)

    if (length(small) == 0L) break

    target <- small[1L]

    donors <- which(counts > min_size)

    if (length(donors) == 0L) {
      stop("Some clusters were too small to determine sensible starting values. Either specify splits manually, or reduce the number of classes.")
    }

    # Find the donor observation closest to
    # the current target-cluster mean.
    target_mean <- mean(yi[splits == target])

    candidates <- which(splits %in% donors)

    selected <- candidates[
      which.min(abs(yi[candidates] - target_mean))
    ]

    splits[selected] <- target
  }

  splits
}

