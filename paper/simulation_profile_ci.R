# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'OpenMx')

lapply(packagenames, library, character.only = TRUE)

outdir <- "paper/ci_profile"
dir.create(outdir)
f <- list.files(outdir)
file.remove(f)
# Read summarydata
summarydata <- readRDS(file = "paper/summarydata.RData")

library(doSNOW)
library(foreach)
cl<-makeCluster((parallel::detectCores()-2L)) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
# summarydata <- summarydata[sample.int(nrow(summarydata), 3, replace = FALSE), ]

max_classes <- max(summarydata$classes)+1L
ci_profile <- foreach(rownum = 1:nrow(summarydata), .packages = packagenames, ,.combine = c, .errorhandling = "remove") %dopar% {
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
  res_mx <- metamixture:::remix(
      yi = df$yi,
      vi = df$vi,
      classes = true_k,
      mu = "varying",
      tau2 = "equal",
      run = TRUE
    )
  sum(sign(res_mx[[true_k]]$output$confidenceIntervals["tau2", c(1,3)] - df$tau2[1])) == 0
}
#Close cluster
stopCluster(cl)
rm(cl)
## stop("End of simulation")
saveRDS(ci_profile, "outdir/ci_profile.RData")

