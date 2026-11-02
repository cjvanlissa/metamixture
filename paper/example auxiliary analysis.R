library(OpenMx)
library(tidySEM)
library(metamixture)
set.seed(7894)
df <- sim_remix_fixed_tau(k = 100, es = 2.5)
res2 <- remix(df$yi, df$vi, classes = 2, tau2 = "equal")
Z <- rnorm(100, mean = df$C)
tmp <- BCH(res2, data = Z)

table_results(tmp)


df_tmp <- df
df_tmp$Z <- rnorm(100, mean = df$C)
library(tidySEM)
pseudo_class(x = res2,
             model = nnet::multinom(class ~ Z, data = data),
             data = df_tmp) -> membership_prediction
saveRDS(res2, "res.RData")
