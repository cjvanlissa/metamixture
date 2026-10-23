#' @importFrom methods hasArg
#' @importFrom stats approx cutree density dist dnorm hclust kmeans rnorm runif
remix_starts <- function(model,
                           splits,
                           ...){
  if(!isTRUE(requireNamespace("OpenMx", quietly = TRUE))) {
    return(NULL)
  }
  stopifnot("mxModel is not a `remix` model." = is_remix(model))
  classes <- num_classes(model)
  if(classes == 1L){
    return(model)
  }
  #df_split <- data.frame(model@matrices$Y$values)
  yi <- model@matrices$Y$values
  if(!hasArg(splits)){
    splits <- try({kmeans(x = yi, centers = classes)$cluster}, silent = TRUE)
    if(inherits(splits, "try-error")){
      message("Could not initialize clusters using K-means, switching to hierarchical clustering.")
      splits <- try({cutree(hclust(dist(yi)), k = classes)}, silent = TRUE)
      if(inherits(splits, "try-error")){
        stop("Could not initialize clusters using hierarchical clustering. Consider using a different clustering method, or imputing missing data.")
      }
    }
  }
  stopifnot("Number of unique values in splits must be identical to the number of latent classes." = length(unique(splits)) == classes)

  tab_split <- table(splits)
  small_cats <- tab_split < 2
  if(any(small_cats)){
    which_small <- which(small_cats)
    atleast <- 2+sum(tab_split[small_cats])
    choose_from <- which(tab_split > atleast)
    if(length(choose_from) == 0) stop("Some clusters were too small to determine sensible starting values. Either specify splits manually, or reduce the number of classes.")
    splits[sample(which(splits %in% choose_from), sum(tab_split[small_cats]))] <- splits[splits %in% names(small_cats)[small_cats]]
  }

  if(!classes == length(unique(splits))){
    stop("Argument 'splits' does not identify a number of groups equal to 'classes'.")
  }
  # Order splits from largest number to smallest
  tab_split <- sort(table(splits), decreasing = TRUE)
  splits <- as.integer(ordered(splits, levels = names(tab_split)))

  for(i in seq_len(classes)){
    res_ma <- dl_estimator(y = model@matrices$Y$values[splits == i], v = model@matrices$V$values[splits == i])
    model@matrices[[paste0("Mu", i)]]$values[1,1] <- res_ma["mu"]
    model@matrices[[paste0("Tau", i)]]$values[1,1] <- res_ma["tau2"]
  }
  p_start <- as.numeric(tab_split) / sum(tab_split)
  model@matrices$Eta$values[1, ] <- p_to_eta(p_start)

  strts <- try(OpenMx::mxRun(model, silent = TRUE, suppressWarnings = TRUE), silent = TRUE)
  if(inherits(strts, "try-error")){
    if(grepl("omxAssignFirstParameters", attr(strts, "condition"), fixed = TRUE)){
      strts <- try(OpenMx::omxAssignFirstParameters(model), silent = TRUE)
      strts <- try(OpenMx::mxRun(strts, silent = TRUE, suppressWarnings = TRUE))
    } else {
      strts <- try(OpenMx::mxTryHard(model, extraTries = 100,
                                          silent = TRUE,
                                          verbose = FALSE,
                                          bestInitsOutput = FALSE))
    }
  }

  if(inherits(strts, "try-error")){
    stop("Could not derive suitable starting values for the ", classes, "-class model.")
  }
  return(strts)
}



dl_estimator <- function (y, v)
{
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
