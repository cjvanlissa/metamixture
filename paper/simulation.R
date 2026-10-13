results_function <- function(df, res){


  ## lazy-blrt
  select_these <- lower.tri(matrix(ncol = ncol(dat), nrow = ncol(dat)))
  dat_stat = cor(dat)[select_these]
  mr_dat <- mr(dat_stat)
  rep_stat <- do.call(rbind, lapply(res, function(thisres) { replicate(reps, {
    sims <- mclust::sim(modelName = thisres$modelName,
                        parameters = thisres$parameters,
                        n = nrow(dat))[, -1, drop = FALSE]
    cr_sim <- cor(sims)[select_these]
    ## srrmr
    sr <- srmr(dat_stat, cr_sim)
    ## mr
    mm <- mr(cr_sim)
    c(mm = mm, srmr = sr)
  })}))

  dif_sr <- rep_stat[seq(from = 4, to = nrow(rep_stat), by = 2), ] - rep_stat[seq(from = 2, to = (nrow(rep_stat)-2), by = 2), ]
  dif_mr <- rep_stat[seq(from = 3, to = nrow(rep_stat), by = 2), ] - rep_stat[seq(from = 1, to = (nrow(rep_stat)-2), by = 2), ]

  ### compile results
  c(
    rowMeans(rep_stat),
    apply(rep_stat, 1, sd),
    as.vector(apply(rep_stat, 1, quantile, probs = c(.025, .975))),
    rowMeans(rep_stat[seq(from = 2, to = nrow(rep_stat), by = 2), ] < .1),
    rowMeans(rep_stat[seq(from = 1, to = nrow(rep_stat), by = 2), ] - mr_dat),
    rowMeans(rep_stat[seq(from = 1, to = nrow(rep_stat), by = 2), ] > mr_dat),
    rowMeans(dif_sr),
    as.vector(apply(dif_sr, 1, quantile, probs = c(.025, .975))),
    rowMeans(dif_sr < 0),
    rowMeans(dif_mr),
    as.vector(apply(dif_mr, 1, quantile, probs = c(.025, .975))),
    rowMeans(dif_mr < 0),
    blrt$p.value,
    sapply(res, AIC),
    sapply(res, BIC),
    ic_ci[seq(1, length(ic_ci), 2)],
    ic_ci[seq(2, length(ic_ci), 2)])
}


# Load packages
packagenames <- c('metamixture', 'OpenMx', 'tidySEM', 'foreach')
lapply(packagenames, library, character.only = TRUE)
num_sims = 10
# Check package versions
# Which packages require correct versions and which versions are they?--------------------------------------------------
versions <- c(
  compareVersion(as.character(packageVersion("pema")), "0.1.0.9000"),
  compareVersion(as.character(packageVersion("rstan")), "2.21.2"))
if(!all(versions > -1)) stop("Using the incorrect version of one or more packages.")

# Prepare seeds
RNGkind("L'Ecuyer-CMRG") # 1. RNG capable of creating multiple streams
set.seed(9148)
seeds <- vector("list", num_sims)
seeds[[1]] <- .Random.seed
for (i in 2:length(seeds)) {
  seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
}
# Create grid with simulation parameters and save it
summarydata <- do.call(rbind, lapply(1:num_sims, function(i){
  # Set seed
  .Random.seed <- seeds[i][[1]]
  data.frame(#Number of studies per dataset
    k = round(exp(runif(1, min = log(20), max = log(100)))),
    classes = sample.int(3, size = 1),
    # vi value
    vi = runif(1, min = .01, max = .1),
    #Effect size
    r2 = runif(1, min = 0, max = .8),
    p_ratio = exp(runif(1, min = log(1), max = log(3)))
  )
}))
summarydata$seed <- seeds
saveRDS(summarydata, file = "summarydata.RData")

sim_data <- function(k, classes, p_ratio, vi, r2){
  ps <- c(p_ratio, rep(1, classes-1L))
  ps <- ps/sum(ps)
  C <- sample(seq_len(classes), size = k, replace = TRUE, prob = ps)
  f_x <- (C-1L)
  sigma_error <- sqrt(var(f_x) * (1 - r2) / r2)
  vi <- runif(k, min = 0, max = vi)
  y <- f_x + rnorm(k, sd = sigma_error + sqrt(vi))
  data.frame(yi = y, vi = vi, C = C)
}



library(doSNOW)
library(foreach)
cl<-makeCluster(30) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

# Run a random sample of test iterations
#summarydata <- summarydata[sample.int(nrow(summarydata), 4, replace = FALSE), ]
f <- list.files("./simulation_2026/results", full.names = TRUE)
file.remove(f)

foreach(rownum = 1:nrow(summarydata), .packages = c("pema", "sn", "rstan", "glmnet"), .final = function(x) NULL) %dopar% {
  # rownum <- sample.int(nrow(summarydata), 1L)
  .Random.seed <- seeds[rownum][[1]]
  # create dataframe for each group
  df <- with(as.list(summarydata[rownum, ]), {
    sim_data(k = k, classes = classes, p_ratio = p_ratio, vi = vi, r2 = r2)
  })
  res <- metamixture::remix(df$yi, df$vi, classes = 1:4)
  fits <- table_fit(res)
  aic <- which.min(fits$AIC)
  bic <- which.min(fits$BIC)
  pmcs <- metamixture:::pmc_remix(res)
  rowSums(sign(pmcs[,c("lb", "ub")]))

  out <- c( c(rownum, reps), results_function(dat, res))
  write.table(x = t(out), file = sprintf("./simulation_2026/results/out_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

#Close cluster
stopCluster(cl)

#Close cluster
stopCluster(cl)
## stop("End of simulation")

save.image(paste0("tab", Sys.Date(), ".RData"))


# Read files --------------------------------------------------------------
library(data.table)

res <- as.data.table(readRDS(list.files(pattern = "summarydata")))
res[, seed := NULL]
res[, reps := NULL]
design_factors <- names(res)[!names(res) %in% c("rownum", "reps")]
res[, (design_factors) := lapply(design_factors, function(thisvar){
  ordered(res[[thisvar]], levels = hyper_parameters[[thisvar]])
})]

max_class <- 6
reps <- 2

# Rename
vars <- c("rownum", "rep",
          paste0("M_", paste0(c("mm", "srmr"), "_", rep(1:max_class, each = reps))), # rowMeans(rep_stat)
          paste0("SD_", paste0(c("mm", "srmr"), "_", rep(1:max_class, each = reps))), # apply(rep_stat, 1, sd)
          paste0(rep(c("lb", "ub"), 6), "_", paste0(rep(c("mm", "srmr"), each = reps), "_", rep(1:max_class, each = 4))), # as.vector(apply(rep_stat, 1, quantile, probs = c(.025, .975)))
          paste0("srmr_.1_", 1:max_class), # rowMeans(rep_stat[seq(from = 2, to = nrow(rep_stat), by = 2), ] < .1)
          paste0("mm_min_dat_", 1:max_class), # rowMeans(rep_stat[seq(from = 1, to = nrow(rep_stat), by = 2), ] - mr_dat)
          paste0("mm_greater_dat_", 1:max_class), # rowMeans(rep_stat[seq(from = 1, to = nrow(rep_stat), by = 2), ] > mr_dat)
          paste0("dif_srmr_", 1:(max_class-1), "_", 2:max_class), # rowMeans(dif_sr)
          paste0("dif_srmr_", rep(c("lb", "ub"), 2), "_", paste0(rep(1:(max_class-1), each = reps), "_", rep(2:max_class, each = reps))), # as.vector(apply(dif_sr, 1, quantile, probs = c(.025, .975)))
          paste0("dif_srmr_neg_", 1:(max_class-1), "_", 2:max_class), # rowMeans(dif_sr < 0)
          paste0("dif_mm_", 1:(max_class-1), "_", 2:max_class), # rowMeans(dif_mr)
          paste0("dif_mm_", rep(c("lb", "ub"), 2), "_", paste0(rep(1:(max_class-1), each = reps), "_", rep(2:max_class, each = reps))), # as.vector(apply(dif_mr, 1, quantile, probs = c(.025, .975)))
          paste0("dif_mm_neg_", 1:(max_class-1), "_", 2:max_class),# rowMeans(dif_mr < 0)
          paste0("blrt_pv_", 1:(max_class-1), "_", 2:max_class), #blrt$p.value,
          paste0("AICS_", 1:max_class), #IC$aics,
          paste0("BICS_", 1:max_class), #IC$bics,
          #"AIC_min", #which.min(IC$aics),
          #"BIC_min", #which.min(IC$bics),
          paste0("AIC_DCI_", 1:(max_class-1), "_", 2:max_class),#ic_ci[,"AIC_0"],
          paste0("BIC_DCI_", 1:(max_class-1), "_", 2:max_class)#ic_ci[,"BIC_0"]
          )

#f <- list.files("results", full.names = TRUE)
#tab <- lapply(f, fread, header = F)
#tab <- rbindlist(tab)
#setorderv(tab, cols = c("V1"), order=1L, na.last=FALSE)
#if(!(tab$V1[1] == 1 & tail(tab$V1, 1) == nrow(res) & length(unique(tab$V1)) == nrow(res))){
#  c(1:nrow(res))[!c(1:nrow(res)) %in% unique(tab$V1)]
#  stop()
#}
colnames(tab) <- vars
#tab[, "rownum" := NULL]
merged <- merge(res, tab, by = "rownum")

fwrite(merged, paste0("sim_results_", Sys.Date(), ".csv"))
saveRDS(merged, paste0("sim_results_", Sys.Date(), ".RData"))
