# Load packages
packagenames <- c('metamixture', 'OpenMx', 'mclust','tidySEM', 'foreach', 'OpenMx')
lapply(packagenames, library, character.only = TRUE)

# Fixed number, random parameters
# num_sims = (parallel::detectCores()-2L)
# RNGkind("L'Ecuyer-CMRG") # 1. RNG capable of creating multiple streams
# set.seed(9148)
# seeds <- vector("list", num_sims)
# seeds[[1]] <- .Random.seed
# for (i in 2:length(seeds)) {
#   seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
# }
# summarydata <- do.call(rbind, lapply(1:num_sims, function(i){
#   # Set seed
#   .Random.seed <- seeds[i][[1]]
#   data.frame(#Number of studies per dataset
#     k = round(exp(runif(1, min = log(200), max = log(300)))),
#     classes = sample.int(3, size = 1),
#     # vi value
#     vi = runif(1, min = .01, max = .1),
#     #Effect size
#     r2 = runif(1, min = 0, max = .8),
#     p_ratio = exp(runif(1, min = log(1), max = log(3)))
#   )
# }))
# summarydata$seed <- seeds


# Factorial design
hyper_parameters <- list(
  reps = 1,
  k = c(20, 50, 100, 300),
  classes = c(2L:3L),
  vi = c(.001, .01, .1),
  r2 = c(.2, .6, .8),
  p_ratio = c(1, 1.5))
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
hyper_parameters <- list(
  reps = 1,
  k = c(20, 50, 100, 300),
  classes = 1L,
  vi = c(.001, .01, .1),
  r2 = NA_real_,
  p_ratio = 1)
summarydata_oneclass <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
summarydata <- rbind(summarydata, summarydata_oneclass)
summarydata$rownum <- 1:nrow(summarydata)

# Prepare seeds
RNGkind("L'Ecuyer-CMRG") # 1. RNG capable of creating multiple streams
set.seed(824)
seeds <- vector("list", nrow(summarydata))
seeds[[1]] <- .Random.seed
for (i in 2:length(seeds)) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}
summarydata$seed <- seeds
saveRDS(summarydata, file = "paper/summarydata.RData")

sim_data <- function(k, classes, p_ratio, vi, r2) {
  ps <- c(p_ratio, rep(1, classes - 1L))
  ps <- ps / sum(ps)
  mu <- seq_len(classes) - 1L
  pop_mean <- sum(ps * mu)
  pop_var <- sum(ps * (mu - pop_mean)^2)

  # Residual heterogeneity, fix for classes = 1 or r2 = 0
  if (classes == 1L) {
    tau2_error <- 1
  } else if (r2 == 0) {
    tau2_error <- 1
  } else {
    tau2_error <- pop_var * (1 - r2) / r2
  }

  C <- sample(seq_len(classes), size = k, replace = TRUE, prob = ps)

  # Sampling variances
  vi_indiv <- runif(k, min = 0, max = vi * tau2_error)

  # Generate observed effects
  y <- mu[C] + rnorm(k, sd = sqrt(tau2_error + vi_indiv))

  data.frame(yi = y, vi = vi_indiv, C = C, tau2 = tau2_error)
}

library(doSNOW)
library(foreach)
cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
#summarydata <- summarydata[sample.int(nrow(summarydata), length(cl), replace = FALSE), ]
f <- list.files("paper/results", full.names = TRUE)
file.remove(f)
write(as.character(Sys.time()), file = "paper/results/time.log", append = TRUE)
max_classes <- max(summarydata$classes)+1L
foreach(rownum = 1:nrow(summarydata), .packages = c("mclust", "metamixture", "OpenMx", "tidySEM"), .errorhandling = "remove", .final = function(x) NULL) %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  # summarydata[rownum, -ncol(summarydata)]

  .Random.seed <- summarydata$seed[rownum][[1]]
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    sim_data(k = k, classes = classes, p_ratio = p_ratio, vi = vi, r2 = r2)
  })
  true_k <- summarydata$classes[rownum]
  true_pi <- summarydata$p_ratio[rownum] / (summarydata$p_ratio[rownum] + true_k - 1L)

  # Mclust
  mc_try <- try({

    res_mc <- lapply(
      seq_len(max_classes),
      mclust::Mclust,
      data = df$yi,
      modelNames = "E"
    )

    # Mclust BIC is reverse-coded
    bic_mc <- which.max(sapply(res_mc, `[[`, "bic"))
    icl_mc <- which.min(sapply(res_mc, metamixture:::icl_mclust))

    pmcs <- metamixture:::pmc_mc(res_mc)
    pmc_mc <- sum(cumprod(
      rowSums(sign(pmcs[seq_len(max_classes - 1L),
                        c("lb", "ub"), drop = FALSE])) == -2
    )) + 1L

    min_mu_mc <- which.min(res_mc[[true_k]]$parameters$mean)

    bias_mu_mc <- res_mc[[true_k]]$parameters$mean[min_mu_mc]
    bias_tau_mc <- (res_mc[[true_k]]$parameters$variance$sigmasq - df$tau2[1])
    bias_pi_mc <- (res_mc[[true_k]]$parameters$pro[min_mu_mc] - true_pi)

    # Bootstrap separately: failure does not invalidate other results
    ses <- try(
      mclust::MclustBootstrap(
        res_mc[[true_k]],
        nboot = 100,
        type = "bs"
      ),
      silent = TRUE
    )

    if (!inherits(ses, "try-error")) {

      ses_summary <- summary(ses, what = "se")

      cov_mu_mc <- res_mc[[true_k]]$parameters$mean[min_mu_mc] +
        c(-1.96, 1.96) * ses_summary$mean[min_mu_mc]
      cov_mu_mc <- sum(sign(cov_mu_mc)) == 0

      cov_tau_mc <- res_mc[[true_k]]$parameters$variance$sigmasq +
        c(-1.96, 1.96) * drop(ses_summary$variance)[1]
      cov_tau_mc <- sum(sign(cov_tau_mc - df$tau2[1])) == 0

      cov_pi_mc <- sum(sign(
        quantile(ses$pro[, min_mu_mc], c(.025, .975), na.rm = TRUE) - true_pi
      )) == 0
      if(true_k == 1L) bias_pi_mc <- cov_pi_mc <- NA

    } else {
      cov_mu_mc <- cov_tau_mc <- cov_pi_mc <- NA
    }
  }, silent = TRUE)

  if (inherits(mc_try, "try-error")) {
    bic_mc <- icl_mc <- pmc_mc <- NA_real_
    cov_mu_mc <- cov_tau_mc <- cov_pi_mc <- NA
    bias_mu_mc <- bias_tau_mc <- bias_pi_mc <- NA
  }


  # OPENMX
  mx_try <- try({
    tstart <- Sys.time()
    res_mx <- lapply(
      seq_len(max_classes),
      metamixture:::remix_mx,
      yi = df$yi,
      vi = df$vi,
      mu = "varying",
      tau2 = "equal",
      run = TRUE
    )
    time_mx <- as.numeric(
      difftime(Sys.time(), tstart, units = "secs")
    )

    conv_mx <- sapply(
      res_mx,
      function(x) metamixture:::mx_convergence(x)$checks_pass
    )
    if(any(!conv_mx)){
      tstart <- Sys.time()
      res_mx[which(!conv_mx)] <- lapply(
        which(!conv_mx), function(i){
          metamixture:::random_starts(res_mx[[i]])
        })
      time_mx <- time_mx + as.numeric(
        difftime(Sys.time(), tstart, units = "secs")
      )
      conv_mx <- sapply(
        res_mx,
        function(x) metamixture:::mx_convergence(x)$checks_pass
      )
    }

    valid_mx <- seq_len(sum(cumprod(conv_mx)))

    if (length(valid_mx) > 0L) {
      aic_mx <- which.min(sapply(res_mx[valid_mx], function(x) {
        tryCatch(
          x@output$Minus2LogLikelihood +
            2 * length(x@output$estimate),
          error = function(e) Inf
        )
      }))

      bic_mx <- which.min(sapply(res_mx[valid_mx], function(x) {
        tryCatch(
          x@output$Minus2LogLikelihood +
            log(x@data$numObs) * length(x@output$estimate),
          error = function(e) Inf
        )
      }))

      icl_mx <- which.min(
        sapply(res_mx[valid_mx], metamixture:::icl_mx)
      )

      if (length(valid_mx) > 1L) {
        pmcs <- metamixture:::pmc_mx(res_mx[valid_mx])
        pmcs <- pmcs[1L:(length(valid_mx)-1L), , drop = FALSE]

        pmc_mx <- sum(cumprod(
          rowSums(sign(pmcs[, c("lb", "ub"), drop = FALSE])) == -2
        )) + 1L
      } else {
        pmc_mx <- 1L
      }

    } else {
      pmc_mx <- aic_mx <- bic_mx <- icl_mx <- NA
    }

    # Parameter estimates and coverage for the true model
    if (true_k %in% valid_mx) {
      m2ll_mx <- res_mx[[true_k]]@output$Minus2LogLikelihood
      min_mu_mx <- which.min(
        sapply(seq_len(true_k), function(i) {
          res_mx[[true_k]][[paste0("Mu", i)]]$values[1, 1]
        })
      )

      cov_mu_mx <- res_mx[[true_k]][[paste0("Mu", min_mu_mx)]]$values[1, 1] +
        c(-1.96, 1.96) *
        res_mx[[true_k]]@output$standardErrors[paste0("mu", min_mu_mx), 1]
      cov_mu_mx <- sum(sign(cov_mu_mx)) == 0

      cov_tau_mx <- res_mx[[true_k]]$Tau1$values[1, 1] +
        c(-1.96, 1.96) *
        res_mx[[true_k]]@output$standardErrors["tau2", 1]
      cov_tau_mx <- sum(sign(cov_tau_mx - df$tau2[1])) == 0

      # Extract mixing proportion
      if(true_k > 1L){
        res_mx_pi <- res_mx[[true_k]]$weights$result[1, min_mu_mx]
        cov_pi_mx <- sum(sign(
          res_mx_pi +
            c(-1.96, 1.96) *
            OpenMx::mxSE("weights", res_mx[[true_k]])[1, min_mu_mx] -
            true_pi
        )) == 0
        bias_pi_mx <- (res_mx_pi - true_pi)
      } else {
        bias_pi_mx <- cov_pi_mx <- NA
      }
      bias_mu_mx <- res_mx[[true_k]][[paste0("Mu", min_mu_mx)]]$values[1,1]
      bias_tau_mx <- (res_mx[[true_k]]$Tau1$values[1,1] - df$tau2[1])
    } else {
      m2ll_mx <- NA
      cov_mu_mx <- cov_tau_mx <- cov_pi_mx <- NA
      bias_mu_mx <- bias_tau_mx <- bias_pi_mx <- NA
    }

  }, silent = TRUE)


  if (inherits(mx_try, "try-error")) {
    time_mx <- m2ll_mx <- NA
    aic_mx <- bic_mx <- icl_mx <- pmc_mx <- NA_real_
    cov_mu_mx <- cov_tau_mx <- cov_pi_mx <- NA
    bias_mu_mx <- bias_tau_mx <- bias_pi_mx <- NA
    valid_mx <- integer(0)
  }


  # C++
  cpp_try <- try({
    tstart <- Sys.time()
    res_cpp <- lapply(
      seq_len(max_classes),
      metamixture:::remix_cpp,
      yi = df$yi,
      vi = df$vi,
      tau2 = "equal"
    )
    time_cpp <- as.numeric(
      difftime(Sys.time(), tstart, units = "secs")
    )

    conv_cpp <- sapply(
      res_cpp,
      function(x) metamixture:::check_convergence_cpp(x)$checks_pass
    )

    valid_cpp <- seq_len(sum(cumprod(conv_cpp)))

    if (length(valid_cpp) > 0L) {

      aic_cpp <- which.min(
        sapply(res_cpp[valid_cpp], `[[`, "AIC")
      )

      bic_cpp <- which.min(
        sapply(res_cpp[valid_cpp], `[[`, "BIC")
      )

      icl_cpp <- which.min(
        sapply(res_cpp[valid_cpp], metamixture:::icl_cpp)
      )

      if (length(valid_cpp) > 1L) {
        pmcs <- metamixture:::pmc_cpp(res_cpp[valid_cpp])
        pmcs <- pmcs[1L:(length(valid_cpp)-1L), , drop = FALSE]
        pmc_cpp <- sum(cumprod(
          rowSums(sign(pmcs[, c("lb", "ub"), drop = FALSE])) == -2
        )) + 1L
      } else {
        pmc_cpp <- 1L
      }
    } else {
      aic_cpp <- bic_cpp <- icl_cpp <- pmc_cpp <- NA
    }

    # Parameter coverage for the true model
    if (true_k %in% valid_cpp) {
      m2ll_cpp <- (-2*res_cpp[[true_k]]$logLik)

      cov_mu_cpp <- res_cpp[[true_k]]$mu[1L] +
        c(-1.96, 1.96) * res_cpp[[true_k]]$se_mu[1L]

      cov_mu_cpp <- sum(sign(cov_mu_cpp)) == 0

      cov_tau_cpp <- res_cpp[[true_k]]$tau2[1] +
        c(-1.96, 1.96) * res_cpp[[true_k]]$se_tau2[1]

      cov_tau_cpp <- sum(sign(cov_tau_cpp - df$tau2[1])) == 0

      cov_pi_cpp <- res_cpp[[true_k]]$pi[1] +
        c(-1.96, 1.96) * res_cpp[[true_k]]$se_pi[1]

      cov_pi_cpp <- sum(sign(cov_pi_cpp - true_pi)) == 0

      bias_mu_cpp <- res_cpp[[true_k]]$mu[1L]
      bias_tau_cpp <- (res_cpp[[true_k]]$tau2[1] - df$tau2[1])
      bias_pi_cpp <- (res_cpp[[true_k]]$pi[1] - true_pi)
      if(true_k == 1L) bias_pi_cpp <- cov_pi_cpp <- NA

    } else{
      m2ll_cpp <- bias_mu_cpp <- bias_tau_cpp <- bias_pi_cpp <- cov_mu_cpp <- cov_tau_cpp <- cov_pi_cpp <- NA
    }

  }, silent = TRUE)



  if (inherits(cpp_try, "try-error")) {
    time_cpp <- NA
    aic_cpp <- bic_cpp <- icl_cpp <- pmc_cpp <- NA_real_
    m2ll_cpp <- bias_mu_cpp <- bias_tau_cpp <- bias_pi_cpp <- cov_mu_cpp <- cov_tau_cpp <- cov_pi_cpp <- NA
    valid_cpp <- integer(0)
  }

  out <- data.frame(
    rownum = summarydata$rownum[rownum],
    true_tau2 = df$tau2[1],
    true_pi = true_pi,

    bic_mc = bic_mc,
    icl_mc = icl_mc,
    pmc_mc = pmc_mc,
    bias_mu_mc = bias_mu_mc,
    bias_tau_mc = bias_tau_mc,
    bias_pi_mc = bias_pi_mc,
    cov_mu_mc = cov_mu_mc,
    cov_tau_mc = cov_tau_mc,
    cov_pi_mc = cov_pi_mc,

    valid_mx = length(valid_mx),
    m2ll_mx = m2ll_mx,
    aic_mx = aic_mx,
    bic_mx = bic_mx,
    icl_mx = icl_mx,
    pmc_mx = pmc_mx,
    time_mx = time_mx,
    bias_mu_mx = bias_mu_mx,
    bias_tau_mx = bias_tau_mx,
    bias_pi_mx = bias_pi_mx,
    cov_mu_mx = cov_mu_mx,
    cov_tau_mx = cov_tau_mx,
    cov_pi_mx = cov_pi_mx,

    valid_cpp = length(valid_cpp),
    m2ll_cpp = m2ll_cpp,
    aic_cpp = aic_cpp,
    bic_cpp = bic_cpp,
    icl_cpp = icl_cpp,
    pmc_cpp = pmc_cpp,
    time_cpp = time_cpp,
    bias_mu_cpp = bias_mu_cpp,
    bias_tau_cpp = bias_tau_cpp,
    bias_pi_cpp = bias_pi_cpp,
    cov_mu_cpp = cov_mu_cpp,
    cov_tau_cpp = cov_tau_cpp,
    cov_pi_cpp = cov_pi_cpp
  )
  saveRDS(out, file = paste0("paper/results/out_", rownum, ".rds"), compress = FALSE)
  #write.table(x = t(out), file = sprintf("paper/results/out_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}
write(as.character(Sys.time()), file = "paper/results/time.log", append = TRUE)
#Close cluster
stopCluster(cl)
rm(cl)
## stop("End of simulation")

#save.image(paste0("tab", Sys.Date(), ".RData"))

# Read files --------------------------------------------------------------
library(data.table)

res <- as.data.table(readRDS(list.files(path = "paper", pattern = "summarydata", full.names = TRUE)))
res[, seed := NULL]

design_factors <- names(res)[!names(res) %in% c("rownum", "reps")]
# res[, (design_factors) := lapply(design_factors, function(thisvar){
#   ordered(res[[thisvar]], levels = hyper_parameters[[thisvar]])
# })]


# Rename

f <- list.files("paper/results", pattern = "\\.rds", full.names = TRUE)
tab <- lapply(f, readRDS)
#names(tab) <- basename(f)
#tab <- lapply(f, fread, header = F)
tab <- as.data.table(do.call(rbind, tab))
#tab <- rbindlist(tab, fill = TRUE)
setorderv(tab, cols = c("rownum"), order=1L, na.last=FALSE)
if(!(tab$rownum[1] == 1 & tail(tab$rownum, 1) == nrow(res) & length(unique(tab$rownum)) == nrow(res))){
  c(1:nrow(res))[!c(1:nrow(res)) %in% unique(tab$rownum)]
  stop()
}
#merged <- merge(res, tab, by = "rownum", all.x = TRUE)
merged <- cbind(res, tab)

fwrite(merged, paste0("sim_results_", Sys.Date(), ".csv"))
saveRDS(merged, paste0("sim_results_", Sys.Date(), ".RData"))
