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
  tau2_ratio = c(2, 4, 6),
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
#res[, seed := NULL]

f <- list.files(outdir, pattern = "\\.rds", full.names = TRUE)
tab <- lapply(f, readRDS)

tab <- as.data.table(do.call(rbind, tab))

setorderv(tab, cols = c("rownum"), order=1L, na.last=FALSE)
setorderv(res, cols = c("rownum"), order=1L, na.last=FALSE)

merged <- cbind(res, tab)

saveRDS(merged, file.path(outdir, paste0("sim_results_", Sys.Date(), ".RData")))

# Analysis

library(data.table)
dat <- readRDS(list.files("paper/metaplus/", pattern = "sim_results", full.names = T))

# Remove duplicated column names
dat <- dat[, !duplicated(names(dat)), with = FALSE]
dat[, conv_mx := is.finite(m2ll_mx)]

# 1. Optimization performance ----------------------------

optimization <- dat[, .(
  n = .N,

  # Number of converged models (1-4)
  mean_conv_mx = mean(conv_mx, na.rm = TRUE),
  mean_conv_mp = mean(conv_mp, na.rm = TRUE),

  # Number of paired likelihood comparisons
  n_ll = sum(is.finite(m2ll_mx) & is.finite(m2ll_mp)),

  # Which method obtained the lower minus2LL?
  better_mx = mean(m2ll_mx < m2ll_mp - 1e-4, na.rm = TRUE),
  better_mp = mean(m2ll_mp < m2ll_mx - 1e-4, na.rm = TRUE),
  equal_ll = mean(abs(m2ll_mx - m2ll_mp) <= 1e-4, na.rm = TRUE),

  # Difference: positive means MX has better fit
  #mean_ll_diff = mean(m2ll_mp - m2ll_mx, na.rm = TRUE),
  #median_ll_diff = median(m2ll_mp - m2ll_mx, na.rm = TRUE),

  # Runtime
  #mean_time_mx = mean(time_mx, na.rm = TRUE),
  #mean_time_mp = mean(time_mp, na.rm = TRUE),
  median_time_mx = median(time_mx, na.rm = TRUE),
  median_time_mp = median(time_mp, na.rm = TRUE)
)]

outdir <- "paper/metaplus/results"
dir.create(outdir)

write.csv(optimization, file.path(outdir, "optimization.csv"), row.names = FALSE)

# Bias

# Reshape bias estimates into long format
bias_long <- melt(
  dat,
  measure.vars = patterns(
    mx = "^bias_(mu|tau|pi)_mx$",
    mp = "^bias_(mu|tau|pi)_mp$"
  ),
  variable.name = "parameter",
  value.name = c("mx", "mp")
)

bias_long[, parameter := c("mu", "tau2", "pi")[parameter]]

bias_table <- bias_long[, .(
  n_mx = sum(!is.na(mx)),
  n_mp = sum(!is.na(mp)),

  mean_bias_mx = mean(mx, na.rm = TRUE),
  mean_bias_mp = mean(mp, na.rm = TRUE),

  mean_abs_bias_mx = mean(abs(mx), na.rm = TRUE),
  mean_abs_bias_mp = mean(abs(mp), na.rm = TRUE),

  rmse_mx = sqrt(mean(mx^2, na.rm = TRUE)),
  rmse_mp = sqrt(mean(mp^2, na.rm = TRUE)),

  # Among paired, nonmissing estimates
  mx_more_accurate = mean(abs(mx) < abs(mp), na.rm = TRUE),
  mp_more_accurate = mean(abs(mp) < abs(mx), na.rm = TRUE)
), by = parameter]

write.csv(bias_table, file.path(outdir, "bias.csv"), row.names = FALSE)

# Coverage
coverage_long <- melt(
  dat,
  measure.vars = patterns(
    mx = "^cov_(mu|tau|pi)_mx$",
    mp = "^cov_(mu|tau|pi)_mp$"
  ),
  variable.name = "parameter",
  value.name = c("mx", "mp")
)

coverage_long[, parameter := c("mu", "tau2", "pi")[parameter]]

coverage_table <- coverage_long[, .(
  n_mx = sum(!is.na(mx)),
  n_mp = sum(!is.na(mp)),

  coverage_mx = if (all(is.na(mx))) NA_real_
  else mean(mx, na.rm = TRUE),

  coverage_mp = if (all(is.na(mp))) NA_real_
  else mean(mp, na.rm = TRUE)

  # Deviation from nominal 95% coverage
  #deviation_mx = if (all(is.na(mx))) NA_real_
  #else mean(mx, na.rm = TRUE) - 0.95,

  #deviation_mp = if (all(is.na(mp))) NA_real_
  #else mean(mp, na.rm = TRUE) - 0.95
), by = parameter]


write.csv(coverage_table, file.path(outdir, "coverage.csv"), row.names = FALSE)

summary_table <- data.table(
  method = c("remix", "metaplus"),

  converged = c(
    mean(dat$conv_mx, na.rm = TRUE),
    mean(dat$conv_mp, na.rm = TRUE)
  ),

  runtime = c(
    median(dat$time_mx, na.rm = TRUE),
    median(dat$time_mp, na.rm = TRUE)
  ),

  MAE_mu = c(
    mean(abs(dat$bias_mu_mx), na.rm = TRUE),
    mean(abs(dat$bias_mu_mp), na.rm = TRUE)
  ),

  MAE_tau2 = c(
    mean(abs(dat$bias_tau_mx), na.rm = TRUE),
    mean(abs(dat$bias_tau_mp), na.rm = TRUE)
  ),

  MAE_pi = c(
    mean(abs(dat$bias_pi_mx), na.rm = TRUE),
    mean(abs(dat$bias_pi_mp), na.rm = TRUE)
  ),


  cov_mu = c(
    mean(dat$cov_mu_mx, na.rm = TRUE),
    mean(dat$cov_mu_mp, na.rm = TRUE)
  ),

  cov_tau2 = c(
    mean(dat$cov_tau_mx, na.rm = TRUE),
    NA_real_
  ),

  cov_pi = c(
    mean(dat$cov_pi_mx, na.rm = TRUE),
    NA_real_
  ),
  t(bias_table[, .SD, .SDcols = grep("more_accurate", names(bias_table))])
)
names(summary_table)[ncol(summary_table)-c(2:0)] <- paste0("Acc. ", c("mu", "tau2", "pi"))

write.csv(summary_table, file.path(outdir, "summary_table.csv"), row.names = FALSE)

# bettermodel <- prop.table(table(merged$m2ll_mx < merged$m2ll_mp))
# cat("We conducted a simulation study with ", nrow(merged), "replications. Simulation conditions were ", deparse(hyper_parameters), ", chosen to align with the assumptions of metaplus.\n\n", file = "paper/results_metaplus.txt")
# cat("In ", worcs::report(bettermodel[2]*100, equals = F), "% of replications, the remix model had a lower minus 2 log likelihood than the metaplus model, indicating that it found a better solution.\n\n", file = "paper/results_metaplus.txt", append = TRUE)
# ps <- c("mu", "tau", "pi")
# merged[, (paste0("diff_", ps)) := lapply(ps, function(p){abs(merged[[paste0("bias_", p, "_mp")]]) - abs(merged[[paste0("bias_", p, "_mx")]])})]
# out <- as.numeric(merged[, lapply(.SD, mean, na.rm = T), .SDcols = grep("diff_", names(merged))])
#
# cat("The difference in bias between the two methods was near-zero for all parameters:", paste0(ps, ": ", formatC(out, digits = 2, format = "f")),  "\n\n", file = "paper/results_metaplus.txt", append = TRUE)
#
# out <- as.numeric(merged[, lapply(.SD, mean, na.rm = T), .SDcols = grep("cov_", names(merged))])
# cat("The methods had near-identical confidence interval coverage for mu, but metaplus did not return confidence intervals for tau2 or pi in any replication. While remix has the advantage of yielding confidence intervals, note that coverage was below nominal levels for all parameters.", paste0(grep("cov_", names(merged), value = T), ": ", formatC(out, digits = 2, format = "f")),  "\n\n", file = "paper/results_metaplus.txt", append = TRUE)
#
#
#
