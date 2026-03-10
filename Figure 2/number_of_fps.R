library(ggplot2)
library(dplyr)

# -------------------------------
# Parameters
# -------------------------------
alpha <- 1
shape <- 4
scale <- 0.8
L_max <- 1.5
beta_max <- 6

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

# Generate fold curve
z_vals <- seq(0.05, 10, length.out = 8000)
curve <- data.frame(
  z = z_vals,
  L = L_of_z(z_vals),
  beta = beta_of_z(z_vals)
) |> filter(L >= 0, L <= L_max, beta >= 0, beta <= beta_max)

# Split at cusp
z_c <- (shape - 1) * scale
lower_branch <- curve |> filter(z <= z_c) |> arrange(L)
upper_branch <- curve |> filter(z >= z_c) |> arrange(desc(L))

# -------------------------------
# Bistable polygon (between folds)
# -------------------------------
bistable_poly <- bind_rows(
  lower_branch,
  upper_branch
)

# -------------------------------
# Top-left triangle to fill region above fold at L=0
# -------------------------------
upper_L0 <- upper_branch[which.min(upper_branch$L), ]
L_top <- approx(upper_branch$beta, upper_branch$L, xout = beta_max)$y

top_left_triangle <- data.frame(
  L = c(0, 0, 0.8221, 0),
  beta = c(beta_max, upper_L0$beta, beta_max, beta_max)
)

# -------------------------------
# Cusp point coordinates
# -------------------------------
L_cusp <- L_of_z(z_c)
beta_cusp <- beta_of_z(z_c)

# -------------------------------
# Plot
# -------------------------------
ggplot() +
  # background: 1 fixed point
  geom_rect(aes(xmin = 0, xmax = L_max, ymin = 0, ymax = beta_max), fill = "grey88") +
  
  # bistable region between folds
  geom_polygon(data = bistable_poly, aes(L, beta), fill = "steelblue", alpha = 0.55) +
  
  # top-left triangle
  geom_polygon(data = top_left_triangle, aes(L, beta), fill = "steelblue", alpha = 0.55) +
  
  # fold curves
  geom_path(data = lower_branch, aes(L, beta), colour = "black", linewidth = 1.2) +
  geom_path(data = upper_branch, aes(L, beta), colour = "black", linewidth = 1.2) +
  
  # labels
  annotate("text", x = 1.5/2, y = 2.5, label = "1 fixed point", size = 5, colour="#413f3f") +
  annotate("text", x = 0.55, y = 5.3, label = "3 fixed points", colour = "#165589", alpha = 1, size = 5) +
  
  # 2 fixed points label with arrow pointing to curve
  #annotate("text", x = 0.5, y = 4.6, label = "2 fixed points", angle = -28, size = 5) +
  #annotate("text", x = 0.95, y = 4.9, label = "2 fixed points", angle = -72, size = 5) +
  
  # cusp point
  geom_point(aes(x = L_cusp, y = beta_cusp), colour = "black", size = 2) +
  annotate("text", x = L_cusp + 0.15, y = beta_cusp, 
           label = "Cusp point", colour = "black", size = 5) +
  
  # axes
  coord_cartesian(xlim = c(0, L_max), ylim = c(2, beta_max)) +
  labs(x = expression("Light level, " * L), y = expression("Social coupling, " * beta)) +
  theme_minimal(base_size = 13)


ggsave("Figure 2/analytic_number_of_fps.pdf", width = 7, height = 7)