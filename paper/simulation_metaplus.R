# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'foreach', 'metaplus')
lapply(packagenames, library, character.only = TRUE)

outdir <- "paper/metaplus"
dir.create(outdir)
f <- list.files(outdir, full.names = TRUE)
file.remove(f)

# Factorial design
hyper_parameters <- list(
  reps = 1:100,
  k = c(20, 50, 100, 300),
  classes = 2L,
  vi = c(.001, .01, .1),
  tau2_ratio = c(1, 2, 6),
  p_ratio = c(1, 1.5)
)
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
summarydata$rownum <- 1:nrow(summarydata)

# Prepare seeds
RNGkind("L'Ecuyer-CMRG") # 1. RNG capable of creating multiple streams
set.seed(6994)
seeds <- vector("list", nrow(summarydata))
seeds[[1]] <- .Random.seed
for (i in 2:length(seeds)) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}
summarydata$seed <- seeds
saveRDS(summarydata, file = file.path(outdir, "summarydata.RData"))

sim_data <- function(k, classes, p_ratio, vi, r2) {

  # Class proportions
  ps <- c(p_ratio, rep(1, classes - 1L))
  ps <- ps / sum(ps)

  # Common mean
  mu <- rep(0, classes)

  # Class-specific heterogeneity variances
  # Class 2 has greater heterogeneity than class 1
  tau2 <- if (classes == 1L) {
    1
  } else {
    c(1, r2)
  }

  # True class memberships
  C <- sample(
    seq_len(classes),
    size = k,
    replace = TRUE,
    prob = ps
  )

  # Study-specific sampling variances
  vi_indiv <- runif(
    k,
    min = 0,
    max = vi * min(tau2)
  )

  # Generate observed effect sizes
  y <- rnorm(
    k,
    mean = mu[C],
    sd = sqrt(vi_indiv + tau2[C])
  )

  data.frame(
    yi = y,
    vi = vi_indiv,
    C = C,
    tau2 = tau2[C]
  )
}
library(doSNOW)
library(foreach)
#cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
cl<-makeCluster(4L) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
summarydata <- summarydata[sample.int(nrow(summarydata), 10*length(cl), replace = FALSE), ]
write(as.character(Sys.time()), file = file.path(outdir, "time.log"), append = TRUE)
max_classes <- max(summarydata$classes)+1L
foreach(rownum = 1:nrow(summarydata), .packages = c("mclust", "metamixture", "OpenMx", "tidySEM", "metaplus"), .errorhandling = "remove", .final = function(x) NULL) %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  # summarydata[rownum, -ncol(summarydata)]

  .Random.seed <- summarydata$seed[rownum][[1]]
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    sim_data(k = k, classes = classes, p_ratio = p_ratio, vi = vi, r2 = tau2_ratio)
  })
  true_k <- summarydata$classes[rownum]
  true_pi <- summarydata$p_ratio[rownum] / (summarydata$p_ratio[rownum] + true_k - 1L)

  # METAPLUS
  time_mp <- m2ll_mp <- NA_real_
  bias_mu_mp <- bias_tau_mp <- bias_pi_mp <- NA_real_
  cov_mu_mp <- cov_tau_mp <- cov_pi_mp <- NA
  conv_mp <- FALSE

  mp_try <- try({

    tstart <- Sys.time()
    fit_mp <- metaplus::metaplus(
      yi = df$yi,
      sei = sqrt(df$vi),
      random = "mixture",
      cores = 1,
      notrials = 20,
      justfit = FALSE
    )
    time_mp <- as.numeric(
      difftime(Sys.time(), tstart, units = "secs")
    )

    m2ll_mp <- -2 * as.numeric(logLik(fit_mp))

    # Extract estimates and confidence intervals
    bias_mu_mp <- fit_mp$results["muhat", "Est."]
    bias_tau_mp <- fit_mp$results["tau2", "Est."]-1
    bias_pi_mp <- (1-fit_mp$results["Outlier prob.", "Est."]) - true_pi

    cov_mu_mp <- (sum(sign(fit_mp$results["muhat", c(2:3)])) == 0)
    cov_tau_mp <- (sum(sign((fit_mp$results["tau2", c(2:3)] - 1))) == 0)
    cov_pi_mp <- (sum(sign((fit_mp$results["Outlier prob.", c(2:3)] - (1-true_pi)))) == 0)

    conv_mp <- is.finite(m2ll_mp)

  }, silent = TRUE)

  if (inherits(mp_try, "try-error")) {
    conv_mp <- FALSE
    time_mp <- m2ll_mp <- NA_real_
  }

  # OPENMX
  mx_try <- try({
    tstart <- Sys.time()
    res_mx <- metamixture::remix(
      yi = df$yi,
      vi = df$vi,
      classes = 2L,
      mu = "equal",
      tau2 = "varying",
      run = TRUE
    )
    time_mx <- as.numeric(
      difftime(Sys.time(), tstart, units = "secs")
    )

    conv_mx <- metamixture:::mx_convergence(res_mx)$checks_pass

    # Parameter estimates and coverage for the true model
    if (conv_mx) {
      m2ll_mx <- res_mx@output$Minus2LogLikelihood
      tau_mx <- diag(res_mx@matrices$S$values)
      small_tau_class <- which.min(tau_mx)

      cov_mu_mx <- res_mx$M$values[1, 1] +
        c(-1.96, 1.96) *
        res_mx@output$standardErrors["mu", 1]
      cov_mu_mx <- sum(sign(cov_mu_mx)) == 0

      cov_tau_mx <- tau_mx[small_tau_class] +
        c(-1.96, 1.96) *
        res_mx@output$standardErrors[paste0("tau2_", small_tau_class), 1]
      cov_tau_mx <- sum(sign(cov_tau_mx - 1)) == 0

      # Extract mixing proportion
      res_mx_pi <- res_mx$weights$result[1, small_tau_class]
      cov_pi_mx <- sum(sign(
        res_mx_pi +
          c(-1.96, 1.96) *
          OpenMx::mxSE("weights", res_mx, silent = TRUE)[1, small_tau_class] -
          true_pi
      )) == 0
      bias_pi_mx <- (res_mx_pi - true_pi)

      bias_mu_mx <- res_mx@matrices$M$values[1, 1]
      bias_tau_mx <- (tau_mx[small_tau_class] - 1)
    } else {
      m2ll_mx <- NA
      cov_mu_mx <- cov_tau_mx <- cov_pi_mx <- NA
      bias_mu_mx <- bias_tau_mx <- bias_pi_mx <- NA
    }

  }, silent = TRUE)


  if (inherits(mx_try, "try-error")) {
    time_mx <- m2ll_mx <- NA
    cov_mu_mx <- cov_tau_mx <- cov_pi_mx <- NA
    bias_mu_mx <- bias_tau_mx <- bias_pi_mx <- NA
    conv_mx <- FALSE
  }


  out <- data.frame(
    rownum = summarydata$rownum[rownum],
    true_tau2 = df$tau2[1],
    true_pi = true_pi,

    conv_mx = conv_mx,
    m2ll_mx = m2ll_mx,
    time_mx = time_mx,
    bias_mu_mx = bias_mu_mx,
    bias_tau_mx = bias_tau_mx,
    bias_pi_mx = bias_pi_mx,
    cov_mu_mx = cov_mu_mx,
    cov_tau_mx = cov_tau_mx,
    cov_pi_mx = cov_pi_mx,
    # metaplus
    conv_mp = conv_mp,
    m2ll_mp = m2ll_mp,
    time_mp = time_mp,
    bias_mu_mp = bias_mu_mp,
    bias_tau_mp = bias_tau_mp,
    bias_pi_mp = bias_pi_mp,
    cov_mu_mp = cov_mu_mp,
    cov_tau_mp = cov_tau_mp,
    cov_pi_mp = cov_pi_mp
  )
  saveRDS(out, file = file.path(outdir, paste0("out_", rownum, ".rds")), compress = FALSE)
}
write(as.character(Sys.time()), file = file.path(outdir, "time.log"), append = TRUE)
#Close cluster
stopCluster(cl)
rm(cl)

# Read files --------------------------------------------------------------
library(data.table)

res <- as.data.table(readRDS(list.files(path = "paper", pattern = "summarydata", full.names = TRUE)))
res[, seed := NULL]

f <- list.files(outdir, pattern = "\\.rds", full.names = TRUE)
tab <- lapply(f, readRDS)

tab <- as.data.table(do.call(rbind, tab))

setorderv(tab, cols = c("rownum"), order=1L, na.last=FALSE)
setorderv(res, cols = c("rownum"), order=1L, na.last=FALSE)
if(!(tab$rownum[1] == 1 & tail(tab$rownum, 1) == nrow(res) & length(unique(tab$rownum)) == nrow(res))){
  c(1:nrow(res))[!c(1:nrow(res)) %in% unique(tab$rownum)]
  stop()
}

merged <- cbind(res, tab)

saveRDS(merged, file.path(outdir, paste0("sim_results_", Sys.Date(), ".RData")))
