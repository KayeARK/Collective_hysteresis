library(ggplot2)
library(dplyr)

# -------------------------------
# Parameters
# -------------------------------
alpha <- 1
L_max <- 2.5
beta_max <- 6

# -------------------------------
# Function to compute fold curve branches
# -------------------------------
compute_fold_branches <- function(shape, scale, alpha = 1, L_max = 1.5, beta_max = 6) {
  
  # Local gamma PDF & CDF functions
  f_gamma_local <- function(z) dgamma(z, shape = shape, scale = scale)
  F_gamma_local <- function(z) pgamma(z, shape = shape, scale = scale)
  
  # Local fold curve functions
  beta_of_z_local <- function(z) 1 / f_gamma_local(z)
  L_of_z_local <- function(z) (z - F_gamma_local(z)/f_gamma_local(z)) / alpha
  
  # Generate fold curve
  z_vals <- seq(0.05, 10, length.out = 8000)
  curve <- data.frame(
    z = z_vals,
    L = L_of_z_local(z_vals),
    beta = beta_of_z_local(z_vals)
  ) |> filter(L >= 0, L <= L_max, beta >= 0, beta <= beta_max)
  
  # Split at cusp
  z_c <- (shape - 1) * scale
  lower_branch <- curve |> filter(z <= z_c) |> arrange(L)
  upper_branch <- curve |> filter(z >= z_c) |> arrange(desc(L))
  
  # Return list with branches and cusp info
  list(
    lower_branch = lower_branch,
    upper_branch = upper_branch,
    z_cusp = z_c,
    L_cusp = L_of_z_local(z_c),
    beta_cusp = beta_of_z_local(z_c)
  )
}

# -------------------------------
# Gamma PDF & CDF
# -------------------------------
f_gamma <- function(z) dgamma(z, shape = shape, scale = scale)
F_gamma <- function(z) pgamma(z, shape = shape, scale = scale)

# -------------------------------
# Fold curve
# -------------------------------
beta_of_z <- function(z) 1 / f_gamma(z)
L_of_z <- function(z) (z - F_gamma(z)/f_gamma(z)) / alpha

# Compute branches using the function
branches1 <- compute_fold_branches(1.5, 0.8, alpha, L_max, beta_max)
lower_branch1 <- branches1$lower_branch
upper_branch1 <- branches1$upper_branch
L_cusp1 <- branches1$L_cusp
beta_cusp1 <- branches1$beta_cusp

# Compute branches using the function
branches2 <- compute_fold_branches(2, 0.8, alpha, L_max, beta_max)
lower_branch2 <- branches2$lower_branch
upper_branch2 <- branches2$upper_branch
L_cusp2 <- branches2$L_cusp
beta_cusp2 <- branches2$beta_cusp

# Compute branches using the function
branches3 <- compute_fold_branches(3, 0.8, alpha, L_max, beta_max)
lower_branch3 <- branches3$lower_branch
upper_branch3 <- branches3$upper_branch
L_cusp3 <- branches3$L_cusp
beta_cusp3 <- branches3$beta_cusp

# Compute branches using the function
branches4 <- compute_fold_branches(4, 0.8, alpha, L_max, beta_max)
lower_branch4 <- branches4$lower_branch
upper_branch4 <- branches4$upper_branch
L_cusp4 <- branches4$L_cusp
beta_cusp4 <- branches4$beta_cusp

# Compute branches using the function
branches5 <- compute_fold_branches(5, 0.8, alpha, L_max, beta_max)
lower_branch5 <- branches5$lower_branch
upper_branch5 <- branches5$upper_branch
L_cusp5 <- branches5$L_cusp
beta_cusp5 <- branches5$beta_cusp

# Compute branches using the function
branches6 <- compute_fold_branches(6, 0.8, alpha, L_max, beta_max)
lower_branch6 <- branches6$lower_branch
upper_branch6 <- branches6$upper_branch
L_cusp6 <- branches6$L_cusp
beta_cusp6 <- branches6$beta_cusp

# -------------------------------
# Plot
# -------------------------------
# Define rainbow color palette (red to violet, darkest to lightest)
rainbow_colors <- c("#8B0000", "#FF4500", "#FFD700", "#32CD32", "#1E90FF", "#9932CC")
shape_values <- c(1.5, 2, 3, 4, 5, 6)

# Create legend data
legend_data <- data.frame(
  x = rep(0.1, 6),
  y = seq(5.5, 3.0, length.out = 6),
  shape = shape_values,
  color = rainbow_colors
)

ggplot() +
  
  # fold curves with different rainbow shades
  geom_path(data = lower_branch1, aes(L, beta), colour = rainbow_colors[1], linewidth = 1.2) +  # Dark red, k=1.5
  geom_path(data = upper_branch1, aes(L, beta), colour = rainbow_colors[1], linewidth = 1.2) +
  geom_path(data = lower_branch2, aes(L, beta), colour = rainbow_colors[2], linewidth = 1.2) +  # Orange red, k=2
  geom_path(data = upper_branch2, aes(L, beta), colour = rainbow_colors[2], linewidth = 1.2) +
  geom_path(data = lower_branch3, aes(L, beta), colour = rainbow_colors[3], linewidth = 1.2) +  # Gold, k=3
  geom_path(data = upper_branch3, aes(L, beta), colour = rainbow_colors[3], linewidth = 1.2) +
  geom_path(data = lower_branch4, aes(L, beta), colour = rainbow_colors[4], linewidth = 1.2) +  # Lime green, k=4
  geom_path(data = upper_branch4, aes(L, beta), colour = rainbow_colors[4], linewidth = 1.2) +
  geom_path(data = lower_branch5, aes(L, beta), colour = rainbow_colors[5], linewidth = 1.2) +  # Dodger blue, k=5
  geom_path(data = upper_branch5, aes(L, beta), colour = rainbow_colors[5], linewidth = 1.2) +
  geom_path(data = lower_branch6, aes(L, beta), colour = rainbow_colors[6], linewidth = 1.2) +  # Dark orchid, k=6
  geom_path(data = upper_branch6, aes(L, beta), colour = rainbow_colors[6], linewidth = 1.2) +
  
  # cusp points with matching rainbow shades
  geom_point(aes(x = L_cusp1, y = beta_cusp1), colour = rainbow_colors[1], size = 2) +
  geom_point(aes(x = L_cusp2, y = beta_cusp2), colour = rainbow_colors[2], size = 2) +
  geom_point(aes(x = L_cusp3, y = beta_cusp3), colour = rainbow_colors[3], size = 2) +
  geom_point(aes(x = L_cusp4, y = beta_cusp4), colour = rainbow_colors[4], size = 2) +
  geom_point(aes(x = L_cusp5, y = beta_cusp5), colour = rainbow_colors[5], size = 2) +
  geom_point(aes(x = L_cusp6, y = beta_cusp6), colour = rainbow_colors[6], size = 2) +
  

  # # labels on cusp points
  # geom_text(aes(x = L_cusp1-0.05, y = beta_cusp1-0.2), 
  #          label = "atop(italic(k)==4.0, italic(l)==0.8)",
  #          colour = "black", size = 3, hjust = 0, parse = TRUE) +
  # geom_text(aes(x = L_cusp2-0.05, y = beta_cusp2-0.2), 
  #          label = "atop(italic(k)==1.5, italic(l)==2.1)", 
  #          colour = "black", size = 3, hjust = 0, parse = TRUE) +
  # geom_text(aes(x = L_cusp3-0.05, y = beta_cusp3-0.2), 
  #          label = "atop(italic(k)==2.0, italic(l)==1.6)", 
  #          colour = "black", size = 3, hjust = 0, parse = TRUE) +
  # geom_text(aes(x = L_cusp4-0.05, y = beta_cusp4-0.2), 
  #          label = "atop(italic(k)==8.0, italic(l)==0.4)", 
  #          colour = "black", size = 3, hjust = 0, parse = TRUE) +
  # geom_text(aes(x = L_cusp5-0.05, y = beta_cusp5-0.2), 
  #          label = "atop(italic(k)==16, italic(l)==0.2)", 
  #          colour = "black", size = 3, hjust = 0, parse = TRUE) +
  # geom_text(aes(x = L_cusp6-0.05, y = beta_cusp6-0.2), 
  #          label = "atop(italic(k)==1.2, italic(l)==2.7)", 
  #          colour = "black", size = 3, hjust = 0, parse = TRUE) +
  # axes
  coord_cartesian(xlim = c(0, L_max), ylim = c(0, beta_max)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression("Light level, " * L), y = expression("Social coupling, " * beta)) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5))


ggsave("Figure 3/analytic_number_of_fps_vary_params_shape.pdf", width = 7, height = 7)

# -------------------------------
# Plot underlying gamma distributions
# -------------------------------

# Create data for plotting gamma PDFs
z_range <- seq(0, 10, length.out = 1000)

# Generate gamma PDF data for each parameter set
gamma_data <- data.frame(
  z = rep(z_range, 6),
  density = c(
    dgamma(z_range, shape = 1.5, scale = 0.8),      # branches1
    dgamma(z_range, shape = 2, scale = 0.8),  # branches2  
    dgamma(z_range, shape = 3, scale = 0.8),      # branches3
    dgamma(z_range, shape = 4, scale = 0.8),      # branches4
    dgamma(z_range, shape = 5, scale = 0.8),    # branches5
    dgamma(z_range, shape = 6, scale = 0.8)   # branches6
  ),
  group = rep(c("k=1.5, θ=0.8", "k=2.0, θ=0.8", "k=3.0, θ=0.8", 
                "k=4.0, θ=0.8", "k=5.0, θ=0.8", "k=6.0, θ=0.8"), each = length(z_range)),
  color_index = rep(c(1, 2, 3, 4, 5, 6), each = length(z_range))
)

# Create data for gamma distribution labels
gamma_labels <- data.frame(
  x = c(max(0, (1.5-1)*0.8) + 0.3,
        max(0, (2-1)*0.8) + 0.3,
        max(0, (3-1)*0.8) + 0.3,
        max(0, (4-1)*0.8) + 0.3,
        max(0, (5-1)*0.8) + 0.3,
        max(0, (6-1)*0.8) + 0.3),
  y = c(dgamma(max(0.01, (1.5-1)*0.8), shape = 1.5, scale = 0.8) + 0.05,
        dgamma(max(0.01, (2-1)*0.8), shape = 2, scale = 0.8) + 0.05,
        dgamma(max(0.01, (3-1)*0.8), shape = 3, scale = 0.8) + 0.05,
        dgamma(max(0.01, (4-1)*0.8), shape = 4, scale = 0.8) + 0.05,
        dgamma(max(0.01, (5-1)*0.8), shape = 5, scale = 0.8) + 0.05,
        dgamma(max(0.01, (6-1)*0.8), shape = 6, scale = 0.8) + 0.05),
  label = c("atop(italic(k)==1.5, theta==0.8)",
            "atop(italic(k)==2.0, theta==0.8)",
            "atop(italic(k)==3.0, theta==0.8)",
            "atop(italic(k)==4.0, theta==0.8)",
            "atop(italic(k)==5.0, theta==0.8)",
            "atop(italic(k)==6.0, theta==0.8)"),
  color = rainbow_colors
)

# Create the gamma distribution plot with legend
p_gamma <- ggplot(gamma_data, aes(x = z, y = density, color = factor(color_index))) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(rainbow_colors, as.character(1:6))) +

  
  labs(x = expression("Threshold, " * theta), 
       y = expression("Probability density, " * Gamma(k, l))) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, max(gamma_data$density) * 1.05))+
  scale_x_continuous(expand = c(0, 0.06)) +
  scale_y_continuous(expand = c(0, 0))

print(p_gamma)
ggsave("Figure 3/gamma_distributions_vary_params_shape.pdf", plot = p_gamma, width = 7, height = 7)