library(ggplot2)

# -------------------------------
# Parameters
# -------------------------------
shape <- 4
scale <- 0.8

# -------------------------------
# Generate gamma distribution data
# -------------------------------
z_range <- seq(0, 10, length.out = 1000)
gamma_data <- data.frame(
  z = z_range,
  density = dgamma(z_range, shape = shape, scale = scale)
)

# Calculate cusp point (vertical line position)
z_cusp <- (shape - 1) * scale+2

# Split data for different shading regions
left_data <- gamma_data[gamma_data$z <= z_cusp, ]
right_data <- gamma_data[gamma_data$z >= z_cusp, ]

# -------------------------------
# Create the plot
# -------------------------------
ggplot() +
  # Left side shading (blue)
  geom_area(data = left_data, aes(x = z, y = density), 
            fill = "#0070C0", alpha = 0.7) +
  
  # Right side shading (red)
  geom_area(data = right_data, aes(x = z, y = density), 
            fill = "#D72727", alpha = 0.7) +
  
  # Main curve line
  geom_line(data = gamma_data, aes(x = z, y = density), 
            color = "black", linewidth = 1) +
  
  # Vertical line at cusp
  geom_vline(xintercept = z_cusp, color = "#FFC000", linewidth = 1.5) +
  
  # Styling
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, max(gamma_data$density) * 1.05)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# Save the plot
ggsave("gamma_distribution_plot.pdf", width = 8, height = 6, bg = "white")