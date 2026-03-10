library(ggplot2)
library(dplyr)


# ------------------------------
# Distribution function
# -------------------------------

make_distribution <- function(dist, params) {

  if (dist == "gamma") {
    shape <- params$shape
    scale <- params$scale

    f <- function(z) dgamma(z, shape = shape, scale = scale)
    F <- function(z) pgamma(z, shape = shape, scale = scale)
  }

  else if (dist == "trunc_gamma") {
    shape <- params$shape
    scale <- params$scale
    Tmax  <- params$Tmax

    Z <- pgamma(Tmax, shape = shape, scale = scale)

    f <- function(z) ifelse(z <= Tmax,
                            dgamma(z, shape, scale) / Z,
                            0)

    F <- function(z) pgamma(pmin(z, Tmax), shape, scale) / Z
  }

  else if (dist == "beta") {
    a <- params$a
    b <- params$b
    Tmax <- params$Tmax

    f <- function(z) ifelse(z <= Tmax,
                            dbeta(z/Tmax, a, b) / Tmax,
                            0)

    F <- function(z) pbeta(pmin(z/Tmax,1), a, b)
  }

  else if (dist == "lognormal") {
    meanlog <- params$meanlog
    sdlog <- params$sdlog

    f <- function(z) dlnorm(z, meanlog, sdlog)
    F <- function(z) plnorm(z, meanlog, sdlog)
  }

  else if (dist == "weibull") {
    shape <- params$shape
    scale <- params$scale

    f <- function(z) dweibull(z, shape, scale)
    F <- function(z) pweibull(z, shape, scale)
  }

  else {
    stop("Unknown distribution")
  }

  list(f = f, F = F)
}


# -------------------------------
# Parameters
# -------------------------------
alpha <- 1
L_max <- 2.5  # Increased to allow curves to extend naturally
beta_max <- 8  # Increased to prevent artificial truncation

# -------------------------------
# Function to compute fold curve branches
# -------------------------------
compute_fold_branches <- function(dist, params,
                                  alpha = 1,
                                  L_max = 1,
                                  beta_max = 6,
                                  z_max = 10) {

  D <- make_distribution(dist, params)
  f <- D$f
  F <- D$F

  beta_of_z <- function(z) 1 / f(z)
  L_of_z    <- function(z) (z - F(z)/f(z)) / alpha

  z_vals <- seq(1e-8, z_max, length.out = 8000)

  # Diagnostic: Check boundary behavior
  z_near_zero <- 1e-6
  z_near_tmax <- params$Tmax - 1e-6
  
  cat("For", paste("a =", params$a, ", b =", params$b), ":\n")
  cat("At z ≈ 0:", "L =", L_of_z(z_near_zero), "\n")
  cat("At z ≈ Tmax:", "L =", L_of_z(z_near_tmax), "\n")
  cat("Minimum L in curve:", min(L_of_z(z_vals), na.rm = TRUE), "\n\n")

  curve <- data.frame(
    z = z_vals,
    L = L_of_z(z_vals),
    beta = beta_of_z(z_vals)
  ) |>
    filter(is.finite(L), is.finite(beta)) |>
    filter(L <= L_max, beta >= 0, beta <= beta_max)  # Removed L >= 0 restriction

  # -------- cusp detection (distribution-agnostic) --------
  dlogf <- function(z) {
    eps <- 1e-5
    (log(f(z+eps)) - log(f(z-eps))) / (2*eps)
  }

  idx <- which.min(abs(dlogf(z_vals)))
  z_c <- z_vals[idx]

  lower_branch <- curve |> filter(z <= z_c) |> arrange(L)
  upper_branch <- curve |> filter(z >= z_c) |> arrange(desc(L))

  list(
    lower_branch = lower_branch,
    upper_branch = upper_branch,
    z_cusp = z_c,
    L_cusp = L_of_z(z_c),
    beta_cusp = beta_of_z(z_c)
  )
}


branches1 <- compute_fold_branches(
  "lognormal",
  list(meanlog = 1, sdlog = 1.5),
  alpha, L_max, beta_max, z_max = 15  # Increased z_max for better coverage
)

branches2 <- compute_fold_branches(
  "lognormal",
  list(meanlog = 1, sdlog = 1),
  alpha, L_max, beta_max, z_max = 15
)

branches3 <- compute_fold_branches(
  "lognormal",
  list(meanlog = 1, sdlog = 0.8),
  alpha, L_max, beta_max, z_max = 15
)

branches4 <- compute_fold_branches(
  "lognormal",
  list(meanlog = 1, sdlog = 0.6),
  alpha, L_max, beta_max, z_max = 15
)

branches5 <- compute_fold_branches(
  "lognormal",
  list(meanlog = 1, sdlog = 0.4),
  alpha, L_max, beta_max, z_max = 15
)

branches6 <- compute_fold_branches(
  "lognormal",
  list(meanlog = 1, sdlog = 0.2),
  alpha, L_max, beta_max, z_max = 15
)


# -------------------------------
# Plot
# -------------------------------
# Define rainbow color palette (red to violet, darkest to lightest)
rainbow_colors <- c("#8B0000", "#FF4500", "#FFD700", "#32CD32", "#1E90FF", "#9932CC")

ggplot() +
  
  # fold curves with different rainbow shades
  geom_path(data = branches1$lower_branch, aes(L, beta), colour = rainbow_colors[1], linewidth = 1.2) +  # Dark red
  geom_path(data = branches1$upper_branch, aes(L, beta), colour = rainbow_colors[1], linewidth = 1.2) +
  geom_path(data = branches2$lower_branch, aes(L, beta), colour = rainbow_colors[2], linewidth = 1.2) +  # Orange red
  geom_path(data = branches2$upper_branch, aes(L, beta), colour = rainbow_colors[2], linewidth = 1.2) +
  geom_path(data = branches3$lower_branch, aes(L, beta), colour = rainbow_colors[3], linewidth = 1.2) +  # Gold
  geom_path(data = branches3$upper_branch, aes(L, beta), colour = rainbow_colors[3], linewidth = 1.2) +
  geom_path(data = branches4$lower_branch, aes(L, beta), colour = rainbow_colors[4], linewidth = 1.2) +  # Lime green
  geom_path(data = branches4$upper_branch, aes(L, beta), colour = rainbow_colors[4], linewidth = 1.2) +
  geom_path(data = branches5$lower_branch, aes(L, beta), colour = rainbow_colors[5], linewidth = 1.2) +  # Dodger blue
  geom_path(data = branches5$upper_branch, aes(L, beta), colour = rainbow_colors[5], linewidth = 1.2) +
  geom_path(data = branches6$lower_branch, aes(L, beta), colour = rainbow_colors[6], linewidth = 1.2) +  # Dark orchid
  geom_path(data = branches6$upper_branch, aes(L, beta), colour = rainbow_colors[6], linewidth = 1.2) +
  
  # cusp points with matching rainbow shades
  geom_point(aes(x = branches1$L_cusp, y = branches1$beta_cusp), colour = rainbow_colors[1], size = 2) +
  geom_point(aes(x = branches2$L_cusp, y = branches2$beta_cusp), colour = rainbow_colors[2], size = 2) +
  geom_point(aes(x = branches3$L_cusp, y = branches3$beta_cusp), colour = rainbow_colors[3], size = 2) +
  geom_point(aes(x = branches4$L_cusp, y = branches4$beta_cusp), colour = rainbow_colors[4], size = 2) +
  geom_point(aes(x = branches5$L_cusp, y = branches5$beta_cusp), colour = rainbow_colors[5], size = 2) +
  geom_point(aes(x = branches6$L_cusp, y = branches6$beta_cusp), colour = rainbow_colors[6], size = 2) +
  

  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 6)) +
  scale_x_continuous(expand = c(0, 0.02)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = expression("Light level, " * L), y = expression("Social coupling, " * beta)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5
  ))


ggsave("Vary distributions/analytic_number_of_fps_vary_params_lognormal.pdf", width = 7, height = 7)

# -------------------------------
# Plot underlying gamma distributions
# -------------------------------

# Create data for plotting gamma PDFs
z_range <- seq(0, 10, length.out = 1000)

# Generate gamma PDF data for each parameter set
gamma_data <- data.frame(
  z = rep(z_range, 6),
  density = c(
    dlnorm(z_range, 1, 1.5),      # branches1
    dlnorm(z_range, 1, 1),  # branches2  
    dlnorm(z_range, 1, 0.8),      # branches3
    dlnorm(z_range, 1, 0.6),      # branches4
    dlnorm(z_range, 1, 0.4),    # branches5
    dlnorm(z_range, 1, 0.2)   # branches6
  ),
  group = rep(c("k=4.0, θ=0.2", "k=4.0, θ=0.4", "k=4.0, θ=0.6", 
                "k=4.0, θ=0.8", "k=4.0, θ=1.0", "k=4.0, θ=1.2"), each = length(z_range)),
  color_index = rep(c(1, 2, 3, 4, 5, 6), each = length(z_range))
)

# Create the gamma distribution plot
p_gamma <- ggplot(gamma_data, aes(x = z, y = density, color = factor(color_index))) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(rainbow_colors, as.character(1:6))) +

  labs(x = expression("Threshold, " * theta), 
       y = expression("Probability density, LogNormal(" * mu * "," * sigma * ")")) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, max(gamma_data$density) * 1.05))+
  scale_x_continuous(expand = c(0, 0.07)) +
  scale_y_continuous(expand = c(0, 0))

print(p_gamma)
ggsave("Vary distributions/lognormal_distributions_vary_params.pdf", plot = p_gamma, width = 7, height = 7)