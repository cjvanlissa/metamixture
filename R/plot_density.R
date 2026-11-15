utils::globalVariables(c(".data"))

plot_density_remix <- function(
    x,
    variables = NULL,
    bw = FALSE,
    conditional = FALSE,
    alpha = .2,
    facet_labels = NULL,
    rawdata = FALSE,
    ...) {

  if(rawdata){
    rugdata <- data.frame(yi = x@matrices$Y$values[, 1, drop = TRUE],
                          vi = x@matrices$V$values[, 1, drop = TRUE])
  } else {
    rugdata <- NULL
  }
  .plot_density_fun_remix(.extract_density_data_remix(x = x, variables = variables), bw = bw, conditional = conditional, alpha = alpha, rawdata = rawdata, data = if (rawdata) rugdata else NULL)
}

.extract_density_data_remix <- function(x, variables = NULL, longform = TRUE) {
  yi <- x@matrices$Y$values[, 1, drop = TRUE]
  vi <- x@matrices$V$values[, 1, drop = TRUE]
  mu <- x@matrices$M$values[1, , drop = TRUE]
  tau2 <- diag(x@matrices$S$values)
  pi <- remix_weights(x)
  n <- length(yi)

  valid <- is.finite(yi) & is.finite(vi)
  yi <- yi[valid]
  vi <- vi[valid]

  if (!length(yi)) {
    stop("No valid observations available.")
  }

  # Plot range
  sds <- sqrt(outer(vi, tau2, "+"))

  x <- seq(
    min(c(yi, mu - 4 * max(sds))),
    max(c(yi, mu + 4 * max(sds))),
    length.out = n
  )

  # Class-specific densities averaged over sampling variances
  densities <- vapply(seq_along(mu), function(k) {
    rowMeans(vapply(vi, function(v) {
      stats::dnorm(
        x,
        mean = mu[k],
        sd = sqrt(tau2[k] + v)
      )
    }, numeric(length(x))))
  }, numeric(length(x)))

  weighted <- sweep(densities, 2, pi, "*")

  out <- data.frame(
    x = rep(x, times = length(mu) + 1L),
    y = c(as.vector(weighted), rowSums(weighted)),
    Class = factor(
      rep(
        c(paste0("class", seq_along(mu)), "Total"),
        each = length(x)
      ),
      levels = c(paste0("class", seq_along(mu)), "Total")
    )
  )

  out
}

.plot_density_fun_remix <- function(
    densities,
    bw = FALSE,
    conditional = FALSE,
    alpha = .2,
    rawdata = FALSE,
    data = NULL) {

  classes <- setdiff(levels(factor(densities$Class)), "Total")

  densities$Class <- factor(
    densities$Class,
    levels = c(classes, "Total")
  )

  if (conditional) {

    densities <- densities[densities$Class != "Total", ]

    density_plot <- ggplot2::ggplot(
      densities,
      ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]],
        fill = .data[["Class"]]
      )
    ) +
      ggplot2::geom_area(position = "fill") +
      ggplot2::labs(x = "Effect size", y = "Class probability")

    if (bw) {
      density_plot <- density_plot +
        ggplot2::scale_fill_grey(start = .2, end = .8)
    } else {
      density_plot <- density_plot +
        ggplot2::scale_fill_manual(
          values = get_palette(length(classes))
        )
    }

  } else {

    densities$alpha <- alpha
    densities$alpha[densities$Class == "Total"] <- 0

    if (bw) {

      density_plot <- ggplot2::ggplot(
        densities,
        ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          linetype = .data[["Class"]]
        )
      ) +
        ggplot2::geom_path() +
        ggplot2::scale_linetype_manual(
          values = c(rep(2, length(classes)), 1)
        )

    } else {

      class_colors <- c(
        get_palette(length(classes)),
        "#000000"
      )
      names(class_colors) <- c(classes, "Total")

      density_plot <- ggplot2::ggplot(
        densities,
        ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          fill = .data[["Class"]],
          colour = .data[["Class"]],
          alpha = .data[["alpha"]]
        )
      ) +
        ggplot2::geom_area(position = "identity") +
        ggplot2::scale_colour_manual(values = class_colors) +
        ggplot2::scale_fill_manual(values = class_colors) +
        ggplot2::scale_alpha_continuous(
          range = c(0, alpha),
          guide = "none"
        )
    }

    density_plot <- density_plot +
      ggplot2::labs(x = "Effect size", y = "Density") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  }

  density_plot <- density_plot + ggplot2::theme_bw()

  if (rawdata) {
    if (is.null(data) || !all(c("yi", "vi") %in% names(data))) {
      stop("Argument 'data' must contain columns 'yi' and 'vi'.")
    }

    rug_data <- data[
      is.finite(data$yi) &
        is.finite(data$vi) &
        data$vi >= 0,
      c("yi", "vi"),
      drop = FALSE
    ]

    if (nrow(rug_data)) {

      rug_data$xmin <- rug_data$yi - 1.96 * sqrt(rug_data$vi)
      rug_data$xmax <- rug_data$yi + 1.96 * sqrt(rug_data$vi)

      # Scale rug positions to the plot's vertical range
      ymax <- if (conditional) 1 else max(densities$y, na.rm = TRUE)

      # Lift intervals above the x-axis and spread them vertically
      set.seed(123)
      rug_data$y <- stats::runif(
        nrow(rug_data),
        min = 0.02 * ymax,
        max = 0.08 * ymax
      )

      density_plot <- density_plot +
        ggplot2::geom_segment(
          data = rug_data,
          ggplot2::aes(
            x = .data[["xmin"]],
            xend = .data[["xmax"]],
            y = .data[["y"]],
            yend = .data[["y"]]
          ),
          inherit.aes = FALSE,
          colour = "grey40",
          linewidth = 1,
          alpha = 0.5
        ) +
        ggplot2::geom_point(
          data = rug_data,
          ggplot2::aes(
            x = .data[["yi"]],
            y = .data[["y"]]
          ),
          inherit.aes = FALSE,
          colour = "black",
          size = 1,
          alpha = 0.7
        )
    }
  }
  density_plot
}
