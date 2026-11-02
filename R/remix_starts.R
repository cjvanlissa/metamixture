remix_starts <- function(model, splits = NULL, ...) {
  stopifnot(
    "mxModel is not a remix model." = is_remix(model)
  )

  classes <- num_classes_remix(model)

  yi <- as.numeric(model$Y$values)
  vi <- as.numeric(model$V$values)

  mu_labels <- model$M$labels[1, ]
  tau_labels <- diag(model$S$labels)

  mu_equal <- all(!is.na(mu_labels)) &&
    length(unique(mu_labels)) == 1L

  tau_equal <- all(!is.na(tau_labels)) &&
    length(unique(tau_labels)) == 1L

  dl_overall <- dl_estimator(y = yi, v = vi)

  # One-class model
  if (classes == 1L) {
    model$M$values[1, 1] <- dl_overall["mu"]
    model$S$values[1, 1] <- dl_overall["tau2"]
    return(model)
  }

  # Choose clustering variable
  if (mu_equal && !tau_equal) {
    cluster_data <- (yi - dl_overall["mu"])^2 /
      (vi + dl_overall["tau2"])
  } else {
    cluster_data <- yi
  }

  # Initialize clusters
  if (is.null(splits)) {

    splits <- tryCatch(
      stats::kmeans(
        x = cluster_data,
        centers = classes,
        nstart = 20
      )$cluster,
      error = function(e) {
        message(
          "K-means initialization failed: ",
          conditionMessage(e),
          ". Switching to hierarchical clustering."
        )
        NULL
      }
    )

    if (is.null(splits)) {

      splits <- tryCatch(
        stats::cutree(
          stats::hclust(
            stats::dist(cluster_data)
          ),
          k = classes
        ),
        error = function(e) {
          stop(
            "Could not initialize clusters: ",
            conditionMessage(e)
          )
        }
      )
    }
  }

  # Validate supplied or estimated clusters
  if (length(splits) != length(yi) ||
      anyNA(splits) ||
      length(unique(splits)) != classes) {
    stop(
      "'splits' must identify exactly ",
      classes,
      " nonempty clusters, with one assignment per observation."
    )
  }

  # Repair small clusters
  splits <- repair_small_clusters(
    splits = splits,
    yi = cluster_data,
    min_size = 2L
  )

  if (length(unique(splits)) != classes) {
    stop(paste0("Could not cluster cases into ", classes, " classes. Supply argument `splits` with cluster assignment."))
  }

  # Order clusters by size
  tab_split <- sort(table(splits), decreasing = TRUE)
  splits <- match(as.character(splits), names(tab_split))

  # Class-specific starting values
  dl_class <- lapply(seq_len(classes), function(i) {
    dl_estimator(y = yi[splits == i], v = vi[splits == i])
  })

  # Means
  if (mu_equal) {
    model$M$values[1, ] <- dl_overall["mu"]
  } else {
    model$M$values[1, ] <- vapply(
      dl_class,
      `[[`,
      numeric(1),
      "mu"
    )
  }
  # Curtail starting values to parameter bounds
  model$M$values[1, ] <- pmax(
    model$M$lbound[1, ],
    pmin(
      model$M$ubound[1, ],
      model$M$values[1, ]
    )
  )
  # Heterogeneity variances
  idx <- cbind(
    seq_len(classes),
    seq_len(classes)
  )
  if (tau_equal) {
    model$S$values[idx] <- dl_overall["tau2"]
  } else {
    model$S$values[idx] <- vapply(
      dl_class,
      `[[`,
      numeric(1),
      "tau2"
    )
  }
  # Curtail starting values to parameter bounds
  model$S$values[idx] <- pmax(
    model$S$lbound[idx],
    pmin(
      model$S$ubound[idx],
      model$S$values[idx]
    )
  )

  # Mixing proportions
  p_start <- as.numeric(tab_split) / length(splits)

  model$Eta$values[1, ] <- p_to_eta(p_start)
  model$Eta$values[1, ] <- pmax(
    model$Eta$lbound[1, ],
    pmin(
      model$Eta$ubound[1, ],
      model$Eta$values[1, ]
    )
  )

  # Synchronize equality-constrained parameters
  OpenMx::omxAssignFirstParameters(model)
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

