## =========================================================
## Bird-singing model with individual gamma-distributed thresholds
## Each bird has its own threshold, sings deterministically when input > threshold
## Loops over multiple values of N and plots them all together
## =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

set.seed(5)

## -----------------------------
## Parameters
## -----------------------------

N_values <- c(20, 50, 100)    # values of N to loop over

T_max    <- 400          # time steps to reach steady state
burn_in  <- 100

alpha <- 1.0
beta  <- 5.0

# Gamma distribution parameters for individual thresholds
k_gamma <- 4.0           # shape parameter for gamma distribution
theta_scale <- 0.8       # scale parameter for gamma distribution

L_upper <- 1.5        # maximum light level for sweep

light_levels <- seq(0, L_upper, length.out = 200)

## -----------------------------
## Kernel definitions
## -----------------------------

G_gaussian <- function(d, sigma = 0.15) {
  exp(-d^2 / (2 * sigma^2))
}

G_exponential <- function(d, alpha = 8) {
  exp(-alpha * d)
}

G_uniform <- function(d) {
  rep(1, length(d))
}

kernels <- list(
  Uniform     = list(fun = G_uniform,     params = list())
  #Gaussian    = list(fun = G_gaussian,    params = list(sigma = 0.15)),
  #Exponential = list(fun = G_exponential, params = list(alpha = 8))
)

## -----------------------------
## Precompute spatial structure (now a function of N)
## -----------------------------

build_K <- function(dist_matrix, kernel, N) {
  G_vals <- do.call(kernel$fun,
                    c(list(as.vector(dist_matrix)), kernel$params))
  G_mat <- matrix(G_vals, nrow = N)
  diag(G_mat) <- 0
  K <- G_mat / rowSums(G_mat)
  K
}

## -----------------------------
## Simulate dynamics — theta_birds passed explicitly
## -----------------------------

simulate_sweep_async <- function(light_seq, K, X_init, theta_birds,
                                 lambda = 1.0, dt = 0.1) {
  N       <- length(X_init)
  X       <- X_init
  steady  <- numeric(length(light_seq))
  N_steps <- ceiling(T_max / dt)

  for (i in seq_along(light_seq)) {
    L <- light_seq[i]

    for (t in 1:N_steps) {
      update_prob <- runif(N)
      to_update   <- which(update_prob < lambda * dt)

      if (length(to_update) > 0) {
        S     <- K %*% X
        input <- alpha * L + beta * S
        X[to_update] <- as.numeric(input[to_update] > theta_birds[to_update])
      }
    }

    steady[i] <- mean(X)
  }

  steady
}

## -----------------------------
## Run simulation for all N values and kernels
## -----------------------------

plot_data <- data.frame()

for (N in N_values) {

  # Per-N spatial structure and thresholds
  theta_birds  <- rgamma(N, shape = k_gamma, scale = theta_scale)
  positions    <- matrix(runif(2 * N), ncol = 2)
  dist_matrix  <- as.matrix(dist(positions))
  K_matrices   <- lapply(kernels, build_K, dist_matrix = dist_matrix, N = N)

  for (k in seq_along(kernels)) {
    K <- K_matrices[[k]]

    ## Upward sweep (start mostly silent)
    X0 <- rbinom(N, 1, 0.05)
    res_up <- simulate_sweep_async(light_levels, K, X0, theta_birds)

    ## Downward sweep (start mostly singing)
    X1 <- rbinom(N, 1, 0.95)
    res_down <- simulate_sweep_async(rev(light_levels), K, X1, theta_birds)
    res_down <- rev(res_down)

    plot_data <- rbind(plot_data, data.frame(
      light_level      = rep(light_levels, 2),
      fraction_singing = c(res_up, res_down),
      kernel           = names(kernels)[k],
      direction        = rep(c("Upward", "Downward"), each = length(light_levels)),
      N                = N
    ))
  }

  cat("Finished N =", N, "\n")
}

plot_data$N <- factor(plot_data$N)

## -----------------------------
## Theoretical analysis with individual thresholds
## -----------------------------

compute_mean_field_gamma <- function(L_vals) {

  all_solutions <- data.frame()

  for (L in L_vals) {

    root_equation <- function(x) {
      if (x < 0 || x > 1) return(Inf)
      F_val <- pgamma(alpha * L + beta * x,
                     shape = k_gamma, scale = theta_scale)
      return(F_val - x)
    }

    x_grid      <- seq(0, 1, length.out = 2000)
    f_vals      <- sapply(x_grid, root_equation)
    sign_changes <- which(diff(sign(f_vals)) != 0)

    roots <- c()
    for (i in sign_changes) {
      try({
        root <- uniroot(root_equation,
                       interval = c(x_grid[i], x_grid[i+1]),
                       tol = 1e-12)$root
        if (root >= 0 && root <= 1) roots <- c(roots, root)
      }, silent = TRUE)
    }

    if (L <= 0.01 && abs(root_equation(0)) < 1e-10) roots <- c(roots, 0)
    roots <- unique(round(roots, 12))

    for (x_root in roots) {
      z      <- alpha * L + beta * x_root
      dF_dz  <- dgamma(z, shape = k_gamma, scale = theta_scale)
      dF_dx  <- beta * dF_dz
      stability <- ifelse(abs(dF_dx) < 1, "stable", "unstable")

      all_solutions <- rbind(all_solutions, data.frame(
        L = L, x = x_root, stability = stability, derivative = dF_dx
      ))
    }
  }

  return(all_solutions)
}

light_levels_2 <- seq(0, L_upper, length.out = 1000)
theory_results <- compute_mean_field_gamma(light_levels_2)

stable_theory   <- theory_results[theory_results$stability == "stable", ]
unstable_theory <- theory_results[theory_results$stability == "unstable", ]

y_threshold_lower <- 0.5
y_threshold_upper <- 0.5

theory_data <- data.frame()

if (nrow(stable_theory) > 0) {
  stable_sorted <- stable_theory[order(stable_theory$L), ]

  lower_points <- stable_sorted[stable_sorted$x < y_threshold_lower, ]
  upper_points <- stable_sorted[stable_sorted$x >= y_threshold_upper, ]

  if (nrow(lower_points) > 0) {
    theory_data <- rbind(theory_data, data.frame(
      light_level = lower_points$L, fraction_singing = lower_points$x,
      branch = "Lower Stable", branch_id = "lower_1", stringsAsFactors = FALSE
    ))
  }

  if (nrow(upper_points) > 0) {
    theory_data <- rbind(theory_data, data.frame(
      light_level = upper_points$L, fraction_singing = upper_points$x,
      branch = "Upper Stable", branch_id = "upper_1", stringsAsFactors = FALSE
    ))
  }
}

if (nrow(unstable_theory) > 0) {
  unstable_sorted <- unstable_theory[order(unstable_theory$L), ]
  theory_data <- rbind(theory_data, data.frame(
    light_level = unstable_sorted$L, fraction_singing = unstable_sorted$x,
    branch = "Unstable", branch_id = "unstable_1", stringsAsFactors = FALSE
  ))
}

## -----------------------------
## Arrows on steepest parts
## -----------------------------

find_steepest_point <- function(x, y, direction) {
  dx <- diff(x)
  dy <- diff(y)
  slopes <- abs(dy / dx)
  max_slope_idx <- which.max(slopes)

  x_arrow <- x[max_slope_idx] + dx[max_slope_idx] / 2
  y_arrow <- y[max_slope_idx] + dy[max_slope_idx] / 2

  arrow_dx <- dx[max_slope_idx] / sqrt(dx[max_slope_idx]^2 + dy[max_slope_idx]^2) * 0.03
  arrow_dy <- dy[max_slope_idx] / sqrt(dx[max_slope_idx]^2 + dy[max_slope_idx]^2) * 0.03

  if (direction == "Downward") {
    arrow_dx <- -arrow_dx
    arrow_dy <- -arrow_dy
  }

  data.frame(x = x_arrow, y = y_arrow,
             xend = x_arrow + arrow_dx, yend = y_arrow + arrow_dy)
}

arrow_data <- data.frame()

for (kernel_name in names(kernels)) {
  for (n_val in levels(plot_data$N)) {

    up_data <- plot_data[plot_data$kernel == kernel_name &
                           plot_data$direction == "Upward" &
                           plot_data$N == n_val, ]
    arrow_up        <- find_steepest_point(up_data$light_level, up_data$fraction_singing, "Upward")
    arrow_up$kernel <- kernel_name
    arrow_up$direction <- "Upward"
    arrow_up$N      <- n_val

    down_data <- plot_data[plot_data$kernel == kernel_name &
                             plot_data$direction == "Downward" &
                             plot_data$N == n_val, ]
    arrow_down        <- find_steepest_point(down_data$light_level, down_data$fraction_singing, "Downward")
    arrow_down$kernel <- kernel_name
    arrow_down$direction <- "Downward"
    arrow_down$N      <- n_val

    arrow_data <- rbind(arrow_data, arrow_up, arrow_down)
  }
}

arrow_data$N <- factor(arrow_data$N)

## -----------------------------
## Colour palette — one colour per N
## -----------------------------

n_colours <- setNames(brewer.pal(max(3, length(N_values)), "Dark2")[seq_along(N_values)],
                      levels(plot_data$N))

## -----------------------------
## Plot
## -----------------------------

p <- ggplot() +

  # Simulation lines — colour by N, linetype by direction
  geom_line(data = plot_data,
            aes(x = light_level, y = fraction_singing,
                color = N, linetype = direction, group = interaction(N, direction)),
            size = 1.2, alpha = 0.9) +

  # Directional arrows
  geom_segment(data = arrow_data,
               aes(x = x, y = y, xend = xend, yend = yend, color = N),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 1.5, alpha = 0.9) +

  # Theoretical stable branches
  geom_line(data = filter(theory_data, branch == "Lower Stable"),
            aes(x = light_level, y = fraction_singing, group = branch_id),
            color = "black", size = 1.4, alpha = 0.7) +

  geom_line(data = filter(theory_data, branch == "Upper Stable",
                          fraction_singing > 1e-6),
            aes(x = light_level, y = fraction_singing, group = branch_id),
            color = "black", size = 1.4, alpha = 0.7) +

  # Theoretical unstable branch
  geom_line(data = filter(theory_data, branch == "Unstable"),
            aes(x = light_level, y = fraction_singing, group = branch_id),
            color = "black", size = 1.4, alpha = 0.7, linetype = "dotted") +

  scale_color_manual(values = n_colours,
                     name = expression(italic(N))) +
  scale_linetype_manual(values = c("Upward" = "solid", "Downward" = "dashed"),
                        name = "Sweep direction") +

  labs(
    x = expression("Light level, " * L),
    y = expression("Steady state proportion of birds singing")
  ) +

  theme_minimal(base_size = 14) +
  theme(
    axis.title  = element_text(size = 16, face = "bold", color = "#2F4F4F"),
    axis.text   = element_text(size = 12, color = "#4F4F4F"),
    legend.title = element_text(size = 12, face = "bold", color = "#2F4F4F"),
    legend.text  = element_text(size = 11, color = "#4F4F4F"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line   = element_line(color = "black", size = 0.5),
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "white", color = "#E0E0E0", size = 0.5)
  ) +

  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1))

print(p)

ggsave(
  filename = "Figure A1/hysteresis_individual_gamma_thresholds_b5_multi_N.pdf",
  plot = p,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  device = "pdf"
)