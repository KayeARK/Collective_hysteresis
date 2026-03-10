library(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(dplyr)

# Load PDE data
load("PDE Model Images/spatial_model_PDE_data.RData")

# ---- Spatial grid ----
Nx <- 100
Ny <- 100
Lx <- 50
Ly <- 50
dx <- Lx / (Nx-1)
dy <- Ly / (Ny-1)
xgrid <- seq(0, Lx, length.out = Nx)
ygrid <- seq(0, Ly, length.out = Ny)

snap_times <- c(0, 2, 5, 10, 20)

# ---- Plot and save each PDE snapshot separately ----
for (k in seq_along(snapshots)) {
  # Create filename
  filename <- sprintf("PDE Model Images/spatial_snapshot_t_%02d.pdf", as.integer(snap_times[k]))
  
  # Create color palette with more contrasting intermediates
  color_palette <- colorRampPalette(c(
    "#6199CE",  # Blue (low activity)
    "#8FA9D3",  # Light blue
    "#BDB9D8",  # Blue-grey
    "#E6C9DD",  # Light purple-grey
    "#F2D7E2",  # Very light pink
    "#F5E5E7",  # Almost white pink
    "#F8D5D0",  # Light salmon
    "#F2B5A4",  # Light coral
    "#E89578",  # Coral
    "#DD754C",  # Orange-red
    "#D46F6C"   # Red (high activity)
  ))(100)
  
  # Create the plot using base R to match the individual plot style
  pdf(filename, width = 10, height = 8)
  par(mar = c(4, 4, 1, 6), bg = "white", family = "serif")
  
  # Set up plot area with exact limits and no padding - force exact user coordinates
  plot.new()
  par(usr = c(0, 50, 0, 50))  # Force exact user coordinates
  
  # Create the heatmap using image() - adjust coordinates to eliminate whitespace
  # Calculate pixel edges to fill exactly 0 to 50
  nx <- ncol(snapshots[[k]])
  ny <- nrow(snapshots[[k]])
  x_centers <- seq(0, 50, length.out = nx)
  y_centers <- seq(0, 50, length.out = ny)
  
  image(x = x_centers, y = y_centers, z = t(snapshots[[k]]), 
        col = color_palette,
        useRaster = TRUE, add = TRUE)
  
  
  # Add a box around the plot
  box(lwd = 2, col = "grey50")
  
  # Add colorbar outside the plot area
  par(xpd = TRUE)
  # Create colorbar as a single image to avoid white lines
  colorbar_matrix <- matrix(seq(0, 1, length.out = 100), nrow = 100, ncol = 1)
  
  # Draw colorbar as an image
  image(x = c(52, 54), y = seq(10, 40, length.out = 100), 
        z = t(colorbar_matrix), col = color_palette, 
        axes = FALSE, useRaster = TRUE, add = TRUE)
  
  # Add colorbar outline and labels
  rect(52, 10, 54, 40, col = NA, border = "black", lwd = 1)
  
  # Add colorbar tick marks and labels
  colorbar_ticks <- c(10, 35/2,25, 65/2,40)  # positions
  colorbar_labels <- c("", "", "", "", "")  # values
  segments(54, colorbar_ticks, 54.5, colorbar_ticks)
  text(55, colorbar_ticks, colorbar_labels, pos = 4, cex = 1)
  
  # Add colorbar title
  text(53, 45, "", srt = 0, pos = 3, cex = 1.2, font = 2)
  
  par(xpd = FALSE)
  
  dev.off()
  
  cat("Saved PDE plot:", filename, "\n")
}

# -------------------------------
# Load and plot Individual data  
# -------------------------------
load("PDE Model Images/spatial_model_individual_data.RData")
individual_snapshots <- snapshots  # Rename to avoid conflicts with PDE data

# Plot individual snapshots
colors <- list(
  inactive = "#6199CE", 
  active = "#D46F6C",
  background = "white"
)

create_spatial_plot <- function(df, filename) {
  pdf(filename, width = 10, height = 8)
  par(mar = c(4, 4, 1, 6), bg = colors$background, family = "serif")
  
  # Set up plot area with exact limits and no padding
  plot.new()
  par(usr = c(0, 50, 0, 50))  # Force exact user coordinates
  
  # Add points
  points(df$x, df$y, 
         col = ifelse(df$X == 1, colors$active, colors$inactive),
         pch = 16, cex = 1.5)
  
  # Add legend outside the plot area
  par(xpd = TRUE)
  legend(x = 55, y = 50, legend = c("Singing", "Silent"), 
         col = c(colors$active, colors$inactive),
         pch = 16, bg = "white", box.col = "grey80", pt.cex = 1.5, cex = 1.2)
  par(xpd = FALSE)
  
  box(lwd = 2, col = "grey50")
  dev.off()
}

# Generate individual plots
for (i in seq_along(individual_snapshots)) {
  time_val <- individual_snapshots[[i]]$time[1]
  filename <- paste0("PDE Model Images/spatial_snapshot_individual_t", time_val, ".pdf")
  create_spatial_plot(individual_snapshots[[i]], filename)
  cat("Saved individual plot:", filename, "\n")
}

