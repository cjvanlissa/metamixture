library(ggplot2)
library(hexSticker)
x <- seq(-4, 5, length.out = 1000)

dens <- data.frame(
  x = x,
  blue = .8*dnorm(x, mean = -0.7, sd = 1),
  orange = dnorm(x, mean = 0.9, sd = 1)
)

# Plot distributions
p <- ggplot(dens, aes(x = x)) +
  geom_area(
    aes(y = blue),
    fill = "#0875C9",
    alpha = 0.45
  ) +
  geom_area(
    aes(y = orange),
    fill = "#FF8500",
    alpha = 0.45
  ) +
  geom_line(
    aes(y = blue),
    color = "#0875C9",
    linewidth = 1.5
  ) +
  geom_line(
    aes(y = orange),
    color = "#FF8500",
    linewidth = 1.5
  ) +
  geom_segment(
    x = -0.7, xend = -0.7,
    y = 0, yend = max(dens$blue),
    color = "#0875C9",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  geom_segment(
    x = 0.9, xend = 0.9,
    y = 0, yend = dnorm(0.9, 0.9, 1),
    color = "#FF8500",
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    color = "#102B49",
    linewidth = 1
  ) +
  coord_cartesian(
    xlim = c(-3.5, 3.7),
    ylim = c(0, 0.46),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    plot.margin = margin(0, 0, 0, 0)
  )

# Generate hex sticker
hexSticker::sticker(
  subplot = p,
  package = "metamixture",

  # Position of the distributions
  s_x = 1,
  s_y = 1.12,
  s_width = 1.65,
  s_height = 1.05,

  # Package name
  p_x = 1,
  p_y = 0.45,
  p_size = 30,
  p_color = "#102B49",
  p_family = "sans",

  # Hexagon
  h_fill = "#FFFFFF",
  h_color = "#102B49",
  h_size = 2,

  # Output
  filename = "metamixture_hex.png",
  dpi = 600
)
