library(data.table)
library(ggplot2)

input <- 'sim_results_2026-10-23.RData'
outdir <- 'paper/output'
ignore <- c("reps", "rownum", "true_tau", "true_pi")
factors <- c('k', 'classes', 'vi', 'r2', 'p_ratio')
methods <- c('bic', 'icl', 'pmc')
pars <- c('mu', 'tau', 'pi')

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "output.txt")
f <- list.files(outdir, full.names = T)
file.remove(f)

dat <- as.data.table(readRDS(input))
setnames(dat, make.unique(names(dat)))
if ('rownum.1' %in% names(dat)) dat[, rownum.1 := NULL]
dat[, (grep("_(mc|cpp)$", names(dat), value = TRUE)) := NULL]
data.table::setnames(dat, grep("_mx$", names(dat), value = TRUE), gsub("_mx$", "", grep("_mx$", names(dat), value = TRUE)))
# Check convergence
conv_tab <- prop.table(table(dat$valid))

cat("In ", worcs::report(sum(dat$valid >= dat$classes)/nrow(dat), equals = F), " percent of replications, all models up to the true model converged.\n", worcs::report(sum(conv_tab[-length(conv_tab)]), equals = F), " percent of replications had one model or more that did not converge.\n", file = outfile, sep = "")
# summary(lm(as.formula(paste0("time ~ ", paste0(fac, collapse = "+"))), data = dat))
#
# cat("In ", worcs::report(sum(dat$valid >= dat$classes)/nrow(dat), equals = F), " percent of replications, all models up to the true model converged.\n", worcs::report(sum(conv_tab[-length(conv_tab)]), equals = F), " percent of replications had one model or more that did not converge.\n", file = outfile)


# Prepare for analysis
set(dat, which(is.na(dat$r2)), "r2", 1.0)

dat[, paste0("bias_abs_", pars) := lapply(.SD, abs), .SDcols = paste0("bias_", pars)]
dat[, paste0("acc_", methods) := lapply(.SD, function(x) x == classes), .SDcols = methods]
dat[, paste0("over_", methods) := lapply(.SD, function(x) x > classes), .SDcols = methods]
dat[, paste0("under_", methods) := lapply(.SD, function(x) x < classes), .SDcols = methods]

dat[, (factors) := lapply(.SD, function(v) ordered(v, levels = unique(v)[order(as.numeric(unique(v)))])), .SDcols = factors]
levels(dat$r2) <- c("0.2", "0.6", "0.8", "1.0")


cols <- setdiff(names(dat), c(ignore, factors))

summary_stats <- rbindlist(lapply(cols, function(col) {
  x <- dat[[col]]
  x_non_na <- x[!is.na(x)]

  data.table(
    outcome = col,
    pct_non_missing = 100 * length(x_non_na) / nrow(dat),
    mean_val = mean(x_non_na),
    sd_val = sd(x_non_na),
    p2_5 = if (length(x_non_na)) quantile(x_non_na, 0.025) else NA_real_,
    p97_5 = if (length(x_non_na)) quantile(x_non_na, 0.975) else NA_real_
  )
}))


cols <- c("true_tau2", "valid", "m2ll", "time",
          "bias_mu", "bias_tau", "bias_pi", "cov_mu", "cov_tau", "cov_pi",
          "bias_abs_mu", "bias_abs_tau", "bias_abs_pi", "acc_bic", "acc_icl",
          "acc_pmc")


sum_by_factor <- rbindlist(lapply(factors, function(fac) {

  rbindlist(lapply(cols, function(col) {

    result <- dat[, {
      x <- get(col)
      x_non_na <- x[!is.na(x)]

      .(
        complete = 100 * length(x_non_na) / .N,
        mean = mean(x_non_na),
        sd = sd(x_non_na),
        p_10 = if (length(x_non_na))
          as.numeric(quantile(x_non_na, 0.1)) else NA_real_,
        p_90 = if (length(x_non_na))
          as.numeric(quantile(x_non_na, 0.9)) else NA_real_
      )
    }, by = fac]

    setnames(result, fac, "level")
    result[, `:=`(outcome = col, factor = fac)]

    result

  }), fill = TRUE)

}), fill = TRUE)

setcolorder(sum_by_factor, c(
  "outcome", "factor", "level",
  "complete", "mean",
  "sd", "p_10", "p_90"
))

setorder(sum_by_factor, outcome, factor, level)

fwrite(summary_stats, file.path(outdir, '01_overall_performance.csv'))
fwrite(sum_by_factor, file.path(outdir, '02_per_condition.csv'))

# Three plots: bias (signed + absolute), coverage, and enumeration
# Prepare plotting data

base_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )

# 1. Bias and absolute bias
# Reshape bias variables into long format
bias_long <- melt(
  dat,
  id.vars = factors,
  measure.vars = patterns(
    bias = "^bias_(mu|tau|pi)$",
    absolute_bias = "^bias_abs_(mu|tau|pi)$"
  ),
  variable.name = "parameter"
)

bias_long[, parameter := factor(
  parameter,
  levels = 1:3,
  labels = c("mu", "tau", "pi")
)]

# Reshape factors into long format
plotdata <- melt(
  bias_long,
  id.vars = c("parameter", "bias", "absolute_bias"),
  measure.vars = factors,
  variable.name = "factor",
  value.name = "level",
  variable.factor = FALSE
)

# Summarize across all other factors
plotdata <- plotdata[
  ,
  .(
    bias = mean(bias, na.rm = TRUE),
    absolute_bias = mean(absolute_bias, na.rm = TRUE)
  ),
  by = .(factor, level, parameter)
]

# Combine bias and absolute bias into one column
plotdata <- melt(
  plotdata,
  id.vars = c("factor", "level", "parameter"),
  measure.vars = c("bias", "absolute_bias"),
  variable.name = "statistic",
  value.name = "value"
)

plotdata[, statistic := factor(
  statistic,
  levels = c("bias", "absolute_bias"),
  labels = c("Bias", "Absolute bias")
)]

# Preserve the original order of factor levels
factor_levels <- lapply(dat[, ..factors], levels)

ordered_levels <- unique(unlist(factor_levels, use.names = FALSE))

plotdata[, level := factor(
  as.character(level),
  levels = ordered_levels,
  ordered = TRUE
)]

plotdata[, factor := factor(
  factor,
  levels = factors
)]

plotdata[, parameter := factor(
  parameter,
  levels = c("mu", "tau", "pi"),
  labels = c("μ", "τ²", "π")
)]

g_bias <- ggplot(
  plotdata,
  aes(
    x = level,
    y = value,
    color = statistic,
    shape = statistic,
    linetype = statistic,
    group = interaction(parameter, statistic)
  )
) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE, size = 2) +
  scale_color_manual(
    values = c(Bias = "grey30", "Absolute bias" = "grey80")) +
  scale_shape_manual(
    values = c(Bias = 16, "Absolute bias" = 17)) +
  facet_grid(parameter ~ factor, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Mean bias",
       color = "Metric", shape = "Metric", linetype = "Parameter") +
  base_theme

ggsave(file.path(outdir, "plot_bias.pdf"),
       g_bias, width = 15, height = 8)


# 2. Coverage
# Coverage variables and simulation factors
coverage_vars <- c("cov_mu", "cov_tau", "cov_pi")
factors <- c("k", "classes", "vi", "r2", "p_ratio")

# Preserve original factor levels and their ordering
factor_levels <- lapply(dat[, ..factors], levels)

# Reshape coverage variables
coverage_long <- melt(
  dat,
  id.vars = factors,
  measure.vars = coverage_vars,
  variable.name = "outcome",
  value.name = "coverage"
)

coverage_long[, outcome := factor(
  outcome,
  levels = coverage_vars,
  labels = c("μ", "τ²", "π")
)]

# Reshape simulation factors
d <- melt(
  coverage_long,
  id.vars = c("outcome", "coverage"),
  measure.vars = factors,
  variable.name = "condition",
  value.name = "level",
  variable.factor = FALSE
)

# Calculate coverage probability, averaging over other factors
d <- d[
  ,
  .(
    value = if (all(is.na(coverage))) NA_real_
    else mean(coverage, na.rm = TRUE)
  ),
  by = .(outcome, condition, level)
]

# Preserve original factor orders independently
d[, factor_level := paste(condition, as.character(level), sep = "__")]

x_levels <- unlist(
  Map(
    function(f, lev) paste(f, lev, sep = "__"),
    names(factor_levels),
    factor_levels
  ),
  use.names = FALSE
)

d[, factor_level := factor(factor_level, levels = x_levels)]
d[, condition := factor(condition, levels = factors)]

# Plot
g_coverage <- ggplot(
  d,
  aes(
    x = factor_level,
    y = value,
    group = outcome,
    linetype = outcome,
    shape = outcome
  )
) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_line(color = "grey30", na.rm = TRUE) +
  geom_point(color = "grey30", na.rm = TRUE, size = 2) +
  facet_grid(
    ~ condition,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) sub("^.*?__", "", x)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Coverage probability",
    linetype = "Parameter",
    shape = "Parameter"
  ) +
  base_theme

ggsave(file.path(outdir, "plot_coverage.pdf"),
       g_coverage, width = 15, height = 8)

# 3. Class enumeration
library(data.table)

factors <- c("k", "classes", "vi", "r2", "p_ratio")

# Preserve original factor levels
factor_levels <- lapply(dat[, ..factors], levels)

# Reshape enumeration outcomes
enum_long <- melt(
  dat,
  id.vars = factors,
  measure.vars = patterns(
    Accuracy = "^acc_(bic|icl|pmc)$",
    Overfitting = "^over_(bic|icl|pmc)$",
    Underfitting = "^under_(bic|icl|pmc)$"
  ),
  variable.name = "criterion"
)

enum_long[, criterion := factor(
  criterion,
  levels = 1:3,
  labels = c("BIC", "ICL", "PMC")
)]

# Reshape simulation factors
d <- melt(
  enum_long,
  id.vars = c("criterion", "Accuracy", "Overfitting", "Underfitting"),
  measure.vars = factors,
  variable.name = "condition",
  value.name = "level",
  variable.factor = FALSE
)

# Calculate proportions, averaging over all other factors
d <- d[
  ,
  lapply(.SD, function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }),
  by = .(criterion, condition, level),
  .SDcols = c("Accuracy", "Overfitting", "Underfitting")
]

# Reshape metrics into long format
d <- melt(
  d,
  id.vars = c("criterion", "condition", "level"),
  measure.vars = c("Accuracy", "Overfitting", "Underfitting"),
  variable.name = "metric",
  value.name = "value"
)

# Preserve metric and criterion ordering
d[, metric := factor(
  metric,
  levels = c("Accuracy", "Overfitting", "Underfitting")
)]

d[, criterion := factor(
  criterion,
  levels = c("BIC", "ICL", "PMC")
)]

d[, condition := factor(
  condition,
  levels = factors
)]

# Preserve original factor-level ordering independently
d[, factor_level := paste(
  condition, as.character(level), sep = "__"
)]

x_levels <- unlist(
  Map(
    function(f, lev) paste(f, lev, sep = "__"),
    names(factor_levels),
    factor_levels
  ),
  use.names = FALSE
)

d[, factor_level := factor(
  factor_level,
  levels = x_levels
)]

# Remove factor prefixes from x-axis labels
x_labels <- function(x) sub("^.*?__", "", x)

# Plot
g_enum <- ggplot(
  d,
  aes(
    x = factor_level,
    y = value,
    color = metric,
    shape = metric,
    linetype = metric,
    group = metric
  )
) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = c(
    "Accuracy" = "grey15",
    "Overfitting" = "grey50",
    "Underfitting" = "grey75"
  )) +
  scale_shape_manual(values = c(
    "Accuracy" = 16,
    "Overfitting" = 17,
    "Underfitting" = 15
  )) +
  scale_x_discrete(labels = x_labels) +
  facet_grid(
    criterion ~ condition,
    scales = "free_x",
    space = "free_x"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Proportion",
    color = "Metric",
    shape = "Metric",
    linetype = "Metric"
  ) +
  base_theme

ggsave(
  file.path(outdir, "plot_enumeration.pdf"),
  g_enum,
  width = 15,
  height = 9
)


# Descriptive partial eta-squared. Treat simulation factors as categorical.
# Sequential sums of squares depend on term order if the design is unbalanced.
fac <- c("k", "classes", "vi", "r2", "p_ratio")

partial_eta <- function(fit) {
  a <- anova(fit)
  ss <- a[["Sum Sq"]]
  names(ss) <- rownames(a)

  residual_ss <- ss["Residuals"]
  terms <- setdiff(names(ss), "Residuals")

  data.table(
    factor = terms,
    partial_eta2 = as.numeric(
      ss[terms] / (ss[terms] + residual_ss)
    )
  )
}

# Select outcomes, excluding identifiers and factors
outcomes <- c("bias_mu", "bias_tau", "bias_pi", "cov_mu", "cov_tau", "cov_pi", "acc_bic",
              "acc_icl", "acc_pmc")


effects <- list()

for (outcome in outcomes) {

  d <- copy(dat[
    !is.na(get(outcome)),
    c(fac, outcome),
    with = FALSE
  ])

  # Convert logical outcomes to numeric
  if (is.logical(d[[outcome]])) {
    set(d, j = outcome, value = as.numeric(d[[outcome]]))
  }

  # Skip nonnumeric outcomes
  if (!is.numeric(d[[outcome]])) next

  # Remove nonfinite observations
  d <- d[is.finite(get(outcome))]

  if (nrow(d) < 10L ||
      uniqueN(d[[outcome]]) < 2L) next

  # Remove unused factor levels
  d[, (fac) := lapply(.SD, droplevels), .SDcols = fac]

  # Skip factors with fewer than two observed levels
  valid_fac <- fac[
    vapply(d[, ..fac], uniqueN, integer(1)) > 1L
  ]

  if (!length(valid_fac)) next

  # Fit additive model
  fit <- try(
    lm(
      reformulate(valid_fac, response = outcome),
      data = d
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) next

  z <- partial_eta(fit)

  z[, `:=`(
    outcome = outcome,
    n = nobs(fit)
  )]

  effects[[length(effects) + 1L]] <- z
}

# Combine results
effect_sizes <- rbindlist(effects, fill = TRUE)

setcolorder(
  effect_sizes,
  c("outcome", "factor", "partial_eta2", "n")
)
effect_sizes <- dcast(
  effect_sizes,
  outcome ~ factor,
  value.var = "partial_eta2"
)

fwrite(
  effect_sizes,
  file.path(outdir, "05_factor_effect_sizes.csv")
)
