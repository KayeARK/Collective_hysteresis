## =========================================================
## Bird-singing model with individual gamma-distributed thresholds
## Each bird has its own threshold, sings deterministically when input > threshold
## =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

set.seed(1)

## -----------------------------
## Parameters
## -----------------------------

N        <- 4000       # number of birds
T_max    <- 400          # time steps to reach steady state
burn_in  <- 100

alpha <- 1.0
beta  <- 3.0

# Gamma distribution parameters for individual thresholds
k_gamma <- 4.0           # shape parameter for gamma distribution
theta_scale <- 0.8       # scale parameter for gamma distribution

L_upper <- 1.5        # maximum light level for sweep

# Generate individual thresholds for each bird
theta_birds <- rgamma(N, shape = k_gamma, scale = theta_scale)

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
## Individual threshold activation
## -----------------------------

# No activation function needed - we use direct threshold comparison

## -----------------------------
## Precompute spatial structure
## -----------------------------

positions <- matrix(runif(2 * N), ncol = 2)
dist_matrix <- as.matrix(dist(positions))

## -----------------------------
## Build kernel matrix
## -----------------------------

build_K <- function(kernel) {
  G_vals <- do.call(kernel$fun,
                    c(list(as.vector(dist_matrix)), kernel$params))
  
  G_mat <- matrix(G_vals, nrow = N)
  
  # Set diagonal to 0 (no self-interaction)
  diag(G_mat) <- 0
  
  # Normalize rows
  K <- G_mat / rowSums(G_mat)
  K
}

K_matrices <- lapply(kernels, build_K)

## -----------------------------
## Simulate dynamics with individual thresholds
## -----------------------------

simulate_light <- function(L, K) {
  
  X <- rbinom(N, 1, 0.1)
  
  for (t in 1:T_max) {
    S <- K %*% X
    # Each bird compares input to its individual threshold
    input <- alpha * L + beta * S
    X <- as.numeric(input > theta_birds)  # Deterministic: sing if input > threshold
  }
  
  mean(X)
}

## -----------------------------
## Hysteresis sweep with individual thresholds
## -----------------------------

simulate_sweep <- function(light_seq, K, X_init) {
  
  X <- X_init
  steady <- numeric(length(light_seq))
  
  for (i in seq_along(light_seq)) {
    L <- light_seq[i]
    
    for (t in 1:T_max) {
      S <- K %*% X
      # Each bird compares input to its individual threshold
      input <- alpha * L + beta * S
      X <- as.numeric(input > theta_birds)  # Deterministic: sing if input > threshold
    }
    
    steady[i] <- mean(X)
  }
  
  steady
}

simulate_sweep_async <- function(light_seq, K, X_init, lambda = 1.0, dt = 0.1) {
  
  X <- X_init
  steady <- numeric(length(light_seq))
  N_steps <- ceiling(T_max / dt)
  
  for (i in seq_along(light_seq)) {
    L <- light_seq[i]
    
    for (t in 1:N_steps) {
      # Each bird has probability lambda * dt to update
      update_prob <- runif(N)
      to_update <- which(update_prob < lambda * dt)
      
      if (length(to_update) > 0) {
        S <- K %*% X  # Input depends on current state of all birds
        input <- alpha * L + beta * S
        
        # Update selected birds asynchronously
        X[to_update] <- as.numeric(input[to_update] > theta_birds[to_update])
      }
    }
    
    steady[i] <- mean(X)
  }
  
  steady
}

## -----------------------------
## Run simulation for all kernels
## -----------------------------

results_up   <- matrix(NA, length(light_levels), length(kernels))
results_down <- matrix(NA, length(light_levels), length(kernels))

colnames(results_up) <- names(kernels)
colnames(results_down) <- names(kernels)

for (k in seq_along(kernels)) {
  
  K <- K_matrices[[k]]
  
  ## Upward sweep (start mostly silent)
  X0 <- rbinom(N, 1, 0.05)
  results_up[, k] <- simulate_sweep_async(light_levels, K, X0)
  
  ## Downward sweep (start mostly singing)
  X1 <- rbinom(N, 1, 0.95)
  results_down[, k] <- simulate_sweep_async(rev(light_levels), K, X1)
  results_down[, k] <- rev(results_down[, k])
}


#save results for later use
save(results_up, results_down, light_levels, file = "Code/gamma_individual_b3_results.RData")
## -----------------------------
## Theoretical analysis with individual thresholds
## -----------------------------

# For individual thresholds, the mean-field equation becomes:
# x = P(alpha*L + beta*x > theta) = P(theta < alpha*L + beta*x)
# This is the gamma CDF: x = F_gamma(alpha*L + beta*x)
# For finite populations: we use beta*x*N/(N-1) to account for self-interaction


# For uniform kernel (mean-field case), we solve this numerically

compute_mean_field_gamma <- function(L_vals) {
  
  all_solutions <- data.frame()
  
  for (L in L_vals) {
    
    # Define the equation x = F(alpha*L + beta*x)
    # We solve F(alpha*L + beta*x) - x = 0
    root_equation <- function(x) {
      if (x < 0 || x > 1) return(Inf)  # x must be in [0,1]
      F_val <- pgamma(alpha * L + beta * x, 
                     shape = k_gamma, scale = theta_scale)
      return(F_val - x)
    }
    
    # Find all roots by searching for sign changes on a fine grid
    x_grid <- seq(0, 1, length.out = 2000)  # Fine grid
    f_vals <- sapply(x_grid, root_equation)
    
    # Find sign changes indicating roots
    sign_changes <- which(diff(sign(f_vals)) != 0)
    
    roots <- c()
    for (i in sign_changes) {
      try({
        root <- uniroot(root_equation, 
                       interval = c(x_grid[i], x_grid[i+1]),
                       tol = 1e-12)$root
        if (root >= 0 && root <= 1) {
          roots <- c(roots, root)
        }
      }, silent = TRUE)
    }
    
    # Check boundaries only if they are genuine roots
    # For L=0, we should have x=0, but for L>0, x=0 is usually spurious
    if (L <= 0.01 && abs(root_equation(0)) < 1e-10) roots <- c(roots, 0)
    # Don't add x=1 boundary unless it's a genuine solution
    # if (abs(root_equation(1)) < 1e-10) roots <- c(roots, 1)
    
    # Remove duplicates
    roots <- unique(round(roots, 12))
    
    # For each root, compute stability from F'(stuff)
    for (x_root in roots) {
      z <- alpha * L + beta * x_root
      
      # Derivative: d/dx[F(alpha*L + beta*x)] = F'(alpha*L + beta*x) * beta
      dF_dz <- dgamma(z, shape = k_gamma, scale = theta_scale)  # F'(z)
      dF_dx <- beta * dF_dz  # d/dx[F(alpha*L + beta*x)]
      
      # Stability: fixed point is stable if |dF/dx| < 1
      stability <- ifelse(abs(dF_dx) < 1, "stable", "unstable")
      
      all_solutions <- rbind(all_solutions, data.frame(
        L = L,
        x = x_root,
        stability = stability,
        derivative = dF_dx
      ))
    }
  }
  
  return(all_solutions)
}

# Compute theoretical curve for uniform kernel
cat("Computing theoretical mean-field solution with finite-size correction...\n")

light_levels_2 <- seq(0, L_upper, length.out = 1000)  # Restrict to biological range [0,1]
theory_results <- compute_mean_field_gamma(light_levels_2)

# Separate stable and unstable branches
stable_theory <- theory_results[theory_results$stability == "stable", ]
unstable_theory <- theory_results[theory_results$stability == "unstable", ]

cat("Total points:", nrow(theory_results), "(", nrow(stable_theory), "stable,", nrow(unstable_theory), "unstable )\n")

# Create theoretical data for plotting - SIMPLIFIED threshold-based branches
theory_data <- data.frame()

# SIMPLE THRESHOLD-BASED BRANCH ASSIGNMENT
# Adjust these thresholds as needed on a case-by-case basis:
# - Increase y_threshold to separate branches more strictly  
# - Decrease y_threshold to capture more points in upper branch
y_threshold_lower <- 1  # Below this = Lower branch
y_threshold_upper <- 1  # Above this = Upper branch

cat("Using y threshold =", y_threshold_lower, "to separate stable branches\n")

if (nrow(stable_theory) > 0) {
  # Sort all stable points
  stable_sorted <- stable_theory[order(stable_theory$L), ]
  
  # Assign branches based on simple y thresholds (fraction_singing)
  lower_points <- stable_sorted[stable_sorted$x < y_threshold_lower, ]
  upper_points <- stable_sorted[stable_sorted$x >= y_threshold_upper, ]
  
  # Add lower branch
  if (nrow(lower_points) > 0) {
    theory_data <- rbind(theory_data, data.frame(
      light_level = lower_points$L,
      fraction_singing = lower_points$x,
      branch = "Lower Stable",
      branch_id = "lower_1",
      stringsAsFactors = FALSE
    ))
  }
  
  # Add upper branch  
  if (nrow(upper_points) > 0) {
    theory_data <- rbind(theory_data, data.frame(
      light_level = upper_points$L,
      fraction_singing = upper_points$x,
      branch = "Upper Stable",
      branch_id = "upper_1",
      stringsAsFactors = FALSE
    ))
  }
}

# Add unstable points as one continuous branch
if (nrow(unstable_theory) > 0) {
  unstable_sorted <- unstable_theory[order(unstable_theory$L), ]
  theory_data <- rbind(theory_data, data.frame(
    light_level = unstable_sorted$L,
    fraction_singing = unstable_sorted$x,
    branch = "Unstable",
    branch_id = "unstable_1",
    stringsAsFactors = FALSE
  ))
}

## -----------------------------
## Plot hysteresis with ggplot2
## -----------------------------

# Prepare simulation data for ggplot
plot_data <- data.frame(
  light_level = rep(light_levels, length(kernels) * 2),
  fraction_singing = c(as.vector(results_up), as.vector(results_down)),
  kernel = rep(rep(names(kernels), each = length(light_levels)), 2),
  direction = rep(c("Upward", "Downward"), each = length(light_levels) * length(kernels))
)

# Function to find steepest point and add arrow data
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
  
  data.frame(
    x = x_arrow,
    y = y_arrow,
    xend = x_arrow + arrow_dx,
    yend = y_arrow + arrow_dy
  )
}


# Create the plot
p <- ggplot() +
  # Add simulation data
  #only add upward sweep for clarity
  geom_line(data = filter(plot_data, direction == "Upward"), 
            aes(x = light_level, y = fraction_singing, 
                color = kernel, linetype = direction),
            size = 1.2, alpha = 1) +
  # Add theoretical curves - use lines with proper grouping to avoid jumps
  geom_line(data = filter(theory_data, branch == "Lower Stable"),
            aes(x = light_level, y = fraction_singing, group = branch_id),
            color = "#2E86AB", size = 1.5, alpha = 1) +
  
  geom_line(data = filter(theory_data, branch == "Upper Stable", fraction_singing > 1e-6),
            aes(x = light_level, y = fraction_singing, group = branch_id),
            color = "#2E86AB", size = 1.5, alpha = 1) +
  
  # Add theoretical curves (unstable branches)
  geom_line(data = filter(theory_data, branch == "Unstable"),
            aes(x = light_level, y = fraction_singing, group = branch_id),
            color = "#2E86AB", size = 1.5, alpha = 1, linetype = "dotted") +

  # Colorblind-friendly palette
  scale_color_manual(values = c(
    "Uniform" = brewer.pal(3, "Dark2")[1]
  )) +
  
  scale_linetype_manual(values = c("Upward" = "solid", "Downward" = "dashed")) +
  
  labs(
    x = expression("Light level, " * L),
    y = expression("Steady state value of proportion of birds singing")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                              color = "#2F4F4F", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, 
                                 color = "#696969", margin = margin(b = 20)),
    axis.title = element_text(size = 16, face = "bold", color = "#2F4F4F"),
    axis.text = element_text(size = 12, color = "#4F4F4F"),
    legend.title = element_text(size = 12, face = "bold", color = "#2F4F4F"),
    legend.text = element_text(size = 11, color = "#4F4F4F"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "white", color = "#E0E0E0", size = 0.5),
    legend.margin = margin(t = 10)
  ) +
  
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1)) +
  guides(
    color = "none",
    linetype = "none"
  )

# Display and save
print(p)

ggsave(
  filename = "hysteresis_individual_gamma_thresholds_beta3.pdf",
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  dpi = 300,
  device = "pdf"
)
