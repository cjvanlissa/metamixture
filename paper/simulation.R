# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'OpenMx')

lapply(packagenames, library, character.only = TRUE)

# Factorial design
hyper_parameters <- list(
  reps = 1:100,
  k = c(20, 50, 100, 300),
  classes = c(2L:3L),
  vi = c(.001, .01, .1),
  es = c(1, 2, 3),
  p_ratio = c(1, 1.5, 2))
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
hyper_parameters <- list(
  reps = 1:100,
  k = c(20, 50, 100, 300),
  classes = 1L,
  vi = c(.001, .01, .1),
  es = 0,
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

sim_data <- function(k, classes, p_ratio, vi, es) {
  if(classes > 1){
    ps <- c(p_ratio, rep(1, classes - 1L))
    ps <- ps / sum(ps)
  } else {
    ps <- 1
  }
  mu <- (seq_len(classes) - 1L) * es
  # Scale everything by fixed tau2:
  tau2_error <- 1
  C <- sample(seq_len(classes), size = k, replace = TRUE, prob = ps)

  vi_indiv <- runif(k, min = 0, max = vi * tau2_error)
  yi <- rnorm(k, mean = mu[C], sd = sqrt(tau2_error + vi_indiv))

  data.frame(
    yi = yi,
    vi = vi_indiv,
    C = C
    #, tau2 = tau2_error
  )
}

library(doSNOW)
library(foreach)
cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
# summarydata <- summarydata[sample.int(nrow(summarydata), 3, replace = FALSE), ]
f <- list.files("paper/results", full.names = TRUE)
file.remove(f)
saveRDS(summarydata, file = "paper/summarydata.RData")
write(as.character(Sys.time()), file = "paper/results/time.log", append = TRUE)
max_classes <- max(summarydata$classes)+1L
foreach(rownum = 1:nrow(summarydata), .packages = packagenames, .errorhandling = "remove", .final = function(x) NULL) %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  # summarydata[rownum, -ncol(summarydata)]

  .Random.seed <- summarydata$seed[rownum][[1]]
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    sim_data(k = k, classes = classes, p_ratio = p_ratio, vi = vi, es = es)
  })
  true_k <- length(unique(df$C))
  true_pi <- sum(df$C==1)/nrow(df)

  # OPENMX
  mx_try <- try({
    tstart <- Sys.time()
    res_mx <- lapply(
      seq_len(max_classes),
      metamixture:::remix,
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
        pmcs <- metamixture:::pmc_mx2(res_mx[valid_mx])
        pmcs <- pmcs[1L:(length(valid_mx)-1L), , drop = FALSE]

        pmc_mx <- sum(cumprod(
          rowSums(sign(pmcs[, c("lb", "ub"), drop = FALSE])) == -2
        )) + 1L

        cv_mx <- metamixture:::cv_enumeration(res_mx[valid_mx])$best_class
      } else {
        pmc_mx <- 1L
        cv_mx <- 1L
      }

    } else {
      cv_mx <-pmc_mx <- aic_mx <- bic_mx <- icl_mx <- NA
    }

    # Parameter estimates and coverage for the true model
    if (true_k %in% valid_mx) {
      m2ll_mx <- res_mx[[true_k]]@output$Minus2LogLikelihood
      mus <- res_mx[[true_k]]@matrices$M$values[1, ]
      min_mu_mx <- which.min(abs(mus))

      cov_mu_mx <- res_mx[[true_k]]@matrices$M$values[1, min_mu_mx] +
        c(-1.96, 1.96) *
        res_mx[[true_k]]@output$standardErrors[paste0("mu", min_mu_mx), 1]
      cov_mu_mx <- sum(sign(cov_mu_mx)) == 0

      cov_tau_mx <- res_mx[[true_k]]@matrices$S$values[1, 1] +
        c(-1.96, 1.96) *
        res_mx[[true_k]]@output$standardErrors["tau2", 1]
      cov_tau_mx <- sum(sign(cov_tau_mx - df$tau2[1])) == 0

      # Extract mixing proportion
      if(true_k > 1L){
        res_mx_pi <- res_mx[[true_k]]$weights$result[1, min_mu_mx]
        cov_pi_mx <- sum(sign(
          res_mx_pi +
            c(-1.96, 1.96) *
            OpenMx::mxSE("weights", res_mx[[true_k]], silent = TRUE)[1, min_mu_mx] -
            true_pi
        )) == 0
        bias_pi_mx <- (res_mx_pi - true_pi)
      } else {
        bias_pi_mx <- cov_pi_mx <- NA
      }
      bias_mu_mx <- res_mx[[true_k]]@matrices$M$values[1, min_mu_mx]
      bias_tau_mx <- (res_mx[[true_k]]@matrices$S$values[1, 1] - 1)
    } else {
      m2ll_mx <- NA
      cov_mu_mx <- cov_tau_mx <- cov_pi_mx <- NA
      bias_mu_mx <- bias_tau_mx <- bias_pi_mx <- NA
    }

  }, silent = TRUE)


  if (inherits(mx_try, "try-error")) {
    time_mx <- m2ll_mx <- NA
    cv_mx <- aic_mx <- bic_mx <- icl_mx <- pmc_mx <- NA_real_
    cov_mu_mx <- cov_tau_mx <- cov_pi_mx <- NA
    bias_mu_mx <- bias_tau_mx <- bias_pi_mx <- NA
    valid_mx <- NA
  }


  out <- data.frame(
    rownum = summarydata$rownum[rownum],
    true_pi = true_pi,

    valid_mx = length(valid_mx),
    m2ll_mx = m2ll_mx,
    aic_mx = aic_mx,
    bic_mx = bic_mx,
    icl_mx = icl_mx,
    pmc_mx = pmc_mx,
    cv_mx = cv_mx,
    time_mx = time_mx,
    bias_mu_mx = bias_mu_mx,
    bias_tau_mx = bias_tau_mx,
    bias_pi_mx = bias_pi_mx,
    cov_mu_mx = cov_mu_mx,
    cov_tau_mx = cov_tau_mx,
    cov_pi_mx = cov_pi_mx

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
#res[, seed := NULL]

f <- list.files("paper/results", pattern = "\\.rds", full.names = TRUE)
tab <- lapply(f, readRDS)
tab <- as.data.table(do.call(rbind, tab))
setorderv(tab, cols = c("rownum"), order=1L, na.last=FALSE)
setorderv(res, cols = c("rownum"), order=1L, na.last=FALSE)
if(!all(res$rownum == tab$rownum)){
  c(1:nrow(res))[!c(1:nrow(res)) %in% unique(tab$rownum)]
  stop()
}
#merged <- merge(res, tab, by = "rownum", all.x = TRUE)
merged <- cbind(res, tab)

saveRDS(merged, paste0("paper/results/sim_results_", Sys.Date(), ".RData"))

