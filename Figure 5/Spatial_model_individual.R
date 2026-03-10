library(ggplot2)

# -------------------------------
# Parameters
# -------------------------------
alpha <- 1.0
L_ext <- 0.5       # external drive
beta  <- 8.0       # social coupling
sigma2 <- 5      # kernel variance = PDE_sigma2 / 2

shape <- 4         # gamma threshold
scale <- 0.8

snap_times <- c(0, 2, 5, 10, 20)

# -------------------------------
# 2D spatial domain
# -------------------------------
Nx <- 50
Ny <- 50
set.seed(42)
points <- data.frame(
  x = runif(Nx*Ny, 0, 50),
  y = runif(Nx*Ny, 0, 50)
)
Ntot <- nrow(points)

# -------------------------------
# Individual thresholds from Gamma
# -------------------------------
theta <- rgamma(Ntot, shape=shape, scale=scale)

# -------------------------------
# Initial condition: singing in [24,26]x[24,26]
# -------------------------------
X <- rep(0, Ntot)
initial_blob <- (points$x >= 24 & points$x <= 26 & points$y >= 24 & points$y <= 26)
X[initial_blob] <- 1
theta[initial_blob] <- 0  # ensures initial blob stays on

# -------------------------------
# Gaussian kernel K(r,r') function
# -------------------------------
gaussian_kernel <- function(r, rp, sigma2=1) {
  dx <- r[1] - rp[1]
  dy <- r[2] - rp[2]
  exp(-(dx^2 + dy^2)/(2*sigma2)) / (2*pi*sigma2)
}

# -------------------------------
# Precompute normalized kernel matrix K_{ij} with zero diagonal
# -------------------------------
Kmat <- matrix(0, Ntot, Ntot)
for (i in 1:Ntot) {
  ri <- as.numeric(points[i, ])
  for (j in 1:Ntot) {
    if (i != j) {  # zero self-interaction
      rj <- as.numeric(points[j, ])
      Kmat[i,j] <- gaussian_kernel(ri, rj, sigma2)
    }
  }
  # Normalize row to sum = 1 if there are neighbors
  if (sum(Kmat[i, ]) > 0) {
    Kmat[i,] <- Kmat[i,] / sum(Kmat[i,])
  }
}

# -------------------------------
# Time loop (asynchronous Markov updates)
# -------------------------------
dt <- 0.1   # timestep
Tend <- 20
Nt <- ceiling(Tend / dt)

rate <- 1
p_update <- 1 - exp(-rate * dt)

snap_times <- c(0, 2, 5, 10, 20)
snapshots <- list()
snap_index <- 1
snapshots[[snap_index]] <- data.frame(points, X=X, time=0)
snap_index <- snap_index + 1

for (t in 1:Nt) {
  
  # Shuffle order of updates for true asynchronous behavior
  update_order <- sample(Ntot)
  
  for (i in update_order) {
    if (runif(1) < p_update) {
      # compute input based on current state of neighbors only
      input_i <- alpha*L_ext + beta * sum(Kmat[i, ] * X)
      X[i] <- as.numeric(input_i >= theta[i])
    }
  }
  
  # Save snapshots at requested times
  t_now <- t * dt
  if (snap_index <= length(snap_times) && t_now >= snap_times[snap_index]) {
    snapshots[[snap_index]] <- data.frame(points, X=X, time=snap_times[snap_index])
    snap_index <- snap_index + 1
  }
}

#save the underlying data
save(snapshots, file = "Figure 5/spatial_model_individual_data.RData")

