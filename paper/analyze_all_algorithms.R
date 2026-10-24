# Mixture simulation analysis: mclust, OpenMx, C++
# Input: dat.RData saved with saveRDS(), or an RDS file supplied below.
# Outputs: CSV tables and PDF figures in outdir.

library(data.table)
library(ggplot2)
library(rpart)

input <- "sim_results_2026-10-23.RData"
outdir <- "paper/output_all"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

raw <- readRDS(input)
dat <- as.data.table(raw)
setnames(dat, make.unique(names(dat)))
dat[, rownum.1 := NULL]
factors <- c("k", "classes", "vi", "r2", "p_ratio")
algs <- c("mc", "mx")
methods <- c("bic", "icl", "pmc")
pars <- c("mu", "tau", "pi")
stopifnot(all(factors %in% names(dat)))

for (v in factors) dat[, (paste0("f_", v)) := factor(get(v))]
fac <- paste0("f_", factors)

# Long-format outcomes. Each row = simulation replication x algorithm x outcome.
records <- list()
add <- function(algorithm, outcome, metric, value, available = TRUE) {
  if (!available) return(invisible(NULL))
  records[[length(records) + 1L]] <<- data.table(
    rownum = dat$rownum, algorithm = algorithm,
    outcome = outcome, metric = metric, value = as.numeric(value)
  )
}

for (a in algs) {
  validcol <- paste0("valid_", a)
  ok <- rep(TRUE, nrow(dat))
  if (a %in% c("mx", "cpp")) {
    if (!validcol %in% names(dat)) stop("Missing ", validcol)
    ok <- !is.na(dat[[validcol]])
  }
  if (validcol %in% names(dat)) {
    add(a, "fit", "success", as.numeric(ok))
  }
  tcol <- paste0("time_", a)
  if (tcol %in% names(dat)) add(a, "runtime", "seconds", dat[[tcol]])

  for (p in pars) {
    bcol <- paste0("bias_", p, "_", a)
    ccol <- paste0("cov_", p, "_", a)
    if (bcol %in% names(dat)) {
      b <- as.numeric(dat[[bcol]])
      b[!ok] <- NA_real_
      add(a, p, "bias", b)
    }
    if (ccol %in% names(dat)) {
      x <- as.numeric(dat[[ccol]])
      x[!ok] <- NA_real_
      add(a, p, "coverage", x)
    }
  }
  for (m in methods) {
    col <- paste0(m, "_", a)
    if (!col %in% names(dat)) next
    pred <- dat[[col]]
    # A missing prediction is an enumeration failure, not a wrong class.
    add(a, m, "enumeration_accuracy", ifelse(is.na(pred), NA_real_, as.numeric(pred == dat$classes)))
    add(a, m, "selected_classes", pred)
  }
}
long <- rbindlist(records)
long <- merge(long, dat[, c("rownum", factors, fac), with = FALSE], by = "rownum", sort = FALSE)

# Overall performance. n_valid is the number of nonmissing observations.
summary_table <- long[metric != "selected_classes", .(
  n_total = .N,
  n_valid = sum(is.finite(value)),
  mean = if (all(!is.finite(value))) NA_real_ else mean(value[is.finite(value)]),
  sd = if (sum(is.finite(value)) < 2) NA_real_ else sd(value[is.finite(value)]),
  Q2.5 = if (all(!is.finite(value))) NA_real_ else unname(quantile(value[is.finite(value)], .025)),
  Q97.5 = if (all(!is.finite(value))) NA_real_ else unname(quantile(value[is.finite(value)], .975))
), by = .(algorithm, outcome, metric)]

fwrite(summary_table, file.path(outdir, "01_overall_performance.csv"))

# Per-condition summaries, including runtime (not plotted).
by_condition <- long[metric != "selected_classes", .(
  n_valid = sum(is.finite(value)),
  mean = if (all(!is.finite(value))) NA_real_ else mean(value[is.finite(value)])
), by = c(factors, "algorithm", "outcome", "metric")]
fwrite(by_condition, file.path(outdir, "02_per_condition.csv"))


# Enumeration: overfitting, underfitting, and accuracy
score_enumeration <- function(truth, pred) {
  observed <- !is.na(truth) & !is.na(pred)

  truth <- truth[observed]
  pred <- pred[observed]

  data.table(
    n = length(truth),
    overfitting = if (length(truth)) mean(pred > truth) else NA_real_,
    underfitting = if (length(truth)) mean(pred < truth) else NA_real_,
    accuracy = if (length(truth)) mean(pred == truth) else NA_real_
  )
}

enum <- rbindlist(lapply(algs, function(a) {
  rbindlist(lapply(methods, function(m) {
    col <- paste0(m, "_", a)
    if (!col %in% names(dat)) return(NULL)

    cbind(
      data.table(algorithm = a, method = m),
      score_enumeration(dat$classes, dat[[col]])
    )
  }))
}))

fwrite(enum, file.path(outdir, "03_enumeration_metrics.csv"))

# Graphs
metrics <- c("bias", "coverage", "enumeration_accuracy")

for (metric_name in metrics) {

  d <- rbindlist(lapply(factors, function(f) {
    x <- paste0("f_", f)

    tmp <- long[
      metric == metric_name,
      .(value = if (all(is.na(value))) NA_real_
        else mean(value, na.rm = TRUE)),
      by = c("algorithm", "outcome", x)
    ]

    setnames(tmp, x, "factor_level")
    tmp[, condition := f]
    tmp
  }), fill = TRUE)

  # Preserve the order of factor levels
  d[, factor_level := factor_level]
  d[, condition := factor(condition, levels = factors)]

  g <- ggplot(
    d,
    aes(
      x = factor_level,
      y = value,
      color = algorithm,
      group = algorithm
    )
  ) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    facet_grid(
      outcome ~ condition,
      scales = "free_x",
      space = "free_x"
    ) +
    labs(
      x = "Simulation condition",
      y = switch(
        metric_name,
        bias = "Mean bias",
        coverage = "Coverage probability",
        enumeration_accuracy = "Class enumeration accuracy"
      ),
      color = "Algorithm",
      title = switch(
        metric_name,
        bias = "Bias by simulation condition",
        coverage = "Coverage by simulation condition",
        enumeration_accuracy = "Class enumeration by simulation condition"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 9),
      legend.position = "bottom"
    )

  if (metric_name != "bias") {
    g <- g + coord_cartesian(ylim = c(0, 1))
  }

  ggsave(
    file.path(outdir, paste0("plot_", metric_name, ".pdf")),
    g,
    width = 15,
    height = 8
  )
}

# Effect sizes: partial eta-squared for factors and algorithm interactions.
# Algorithm is within-replication; use rownum as a blocking term.
# The resulting ANOVA is descriptive, particularly for binary coverage/accuracy.
partial_eta <- function(fit) {
  a <- anova(fit)
  ss <- a[["Sum Sq"]]
  names(ss) <- rownames(a)
  err <- ss["Residuals"]
  z <- ss[names(ss) != "Residuals"]
  data.table(term = names(z), partial_eta2 = as.numeric(z / (z + err)))
}

# Use replicate-level outcomes. Different algorithms are measured on the same dataset.
effect_data <- long[metric %in% c("bias", "coverage", "enumeration_accuracy") &
                      is.finite(value)]
effect_data[, algorithm := factor(algorithm)]
all_effects <- list()
for (ii in seq_len(nrow(unique(effect_data[, .(outcome, metric)])))) {
  key <- unique(effect_data[, .(outcome, metric)])[ii]
  d <- effect_data[outcome == key$outcome & metric == key$metric]
  if (nrow(d) < 10 || uniqueN(d$algorithm) < 2) next
  # This is a factorial fixed-effects analysis; use algorithm interactions.
  rhs <- paste0("algorithm * (", paste(fac, collapse = " + "), ")")
  fit <- try(lm(as.formula(paste("value ~", rhs)), data = d), silent = TRUE)
  if (inherits(fit, "try-error")) next
  z <- partial_eta(fit)
  z[, `:=`(outcome = key$outcome, metric = key$metric)]
  all_effects[[length(all_effects) + 1L]] <- z
}
effects <- rbindlist(all_effects, fill = TRUE)
fwrite(effects, file.path(outdir, "05_factor_effect_sizes.csv"))

# Paired algorithm comparisons for each outcome and metric.
# Differences use only simulation replicates with both algorithms observed.
comparisons <- list()
for (ii in seq_len(nrow(unique(effect_data[, .(outcome, metric)])))) {
  key <- unique(effect_data[, .(outcome, metric)])[ii]
  d <- effect_data[outcome == key$outcome & metric == key$metric,
                   .(rownum, algorithm, value)]
  wide <- dcast(d, rownum ~ algorithm, value.var = "value")
  pairs <- combn(intersect(algs, names(wide)), 2, simplify = FALSE)
  for (pair in pairs) {
    delta <- wide[[pair[1]]] - wide[[pair[2]]]
    delta <- delta[is.finite(delta)]
    comparisons[[length(comparisons) + 1L]] <- data.table(
      outcome = key$outcome, metric = key$metric,
      algorithm_1 = pair[1], algorithm_2 = pair[2],
      n = length(delta), mean_difference = if (length(delta)) mean(delta) else NA_real_,
      se_difference = if (length(delta) > 1) sd(delta) / sqrt(length(delta)) else NA_real_
    )
  }
}
fwrite(rbindlist(comparisons), file.path(outdir, "06_paired_algorithm_differences.csv"))

# Which algorithm is best per condition? Define a single interpretable objective:
# maximize enumeration accuracy for each method separately.
# Ties are explicitly recorded rather than resolved by column order.
best <- copy(by_condition[metric == "enumeration_accuracy", c(factors, "algorithm", "outcome", "mean", "n_valid"), with = FALSE])
best <- best[is.finite(mean) & n_valid > 0]
best[, best_score := max(mean), by = c(factors, "outcome")]
best[, is_best := abs(mean - best_score) < 1e-12]
best_labels <- best[is_best == TRUE,
                    .(best_algorithm = paste(sort(algorithm), collapse = "+"),
                      best_accuracy = first(best_score)), by = c(factors, "outcome")]
fwrite(best_labels, file.path(outdir, "07_best_enumeration_by_condition.csv"))

# Decision trees for each enumeration criterion (only if enough distinct conditions).
for (m in methods) {
  d <- copy(best_labels[outcome == m])
  if (nrow(d) < 10 || uniqueN(d$best_algorithm) < 2) next
  for (f in factors) d[, (f) := factor(get(f))]
  d[, best_algorithm := factor(best_algorithm)]
  fit <- rpart(as.formula(paste("best_algorithm ~", paste(factors, collapse = "+"))),
               data = d, method = "class", control = rpart.control(minbucket = 3, cp = .005))
  pdf(file.path(outdir, paste0("tree_best_", m, ".pdf")), width = 10, height = 7)
  plot(fit, uniform = TRUE, margin = .12)
  text(fit, use.n = TRUE, cex = .8)
  dev.off()
}

# Algorithm-specific fit statuses, without assuming which codes indicate success.
statuses <- rbindlist(lapply("mx", function(a) {
  col <- paste0("valid_", a)
  if (!col %in% names(dat)) return(NULL)
  z <- dat[, .N, by = col]
  setnames(z, col, "status")
  z[, algorithm := a]
  z
}))
fwrite(statuses, file.path(outdir, "08_fit_status_counts.csv"))

cat("Analysis written to:", normalizePath(outdir), "\n")
