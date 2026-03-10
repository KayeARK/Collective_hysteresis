# ==============================
# 2D Reaction–Diffusion simulation
# Gamma threshold distribution
# ==============================

# ---- Parameters ----
alpha  <- 1.0
beta   <- 8
L_ext  <- 0.5
sigma2 <- 5
d      <- 2               # 2D

shape  <- 4
scale  <- 0.8

# diffusion prefactor
D0 <- beta * sigma2 / (2*d)

# ---- Spatial grid ----
Nx <- 100
Ny <- 100
Lx <- 50
Ly <- 50
dx <- Lx / (Nx-1)
dy <- Ly / (Ny-1)
xgrid <- seq(0, Lx, length.out = Nx)
ygrid <- seq(0, Ly, length.out = Ny)

# ---- Time grid ----
dt <- 0.01
Tend <- 20
Nt <- round(Tend/dt)

# ---- Initial condition ----
x <- matrix(0, nrow=Nx, ncol=Ny)
# small activated patch in the middle
x[48:52, 48:52] <- 1

# ---- Helper: 2D Laplacian with Neumann BC ----
laplacian2D <- function(u, dx, dy) {
  lap <- matrix(0, nrow=nrow(u), ncol=ncol(u))
  
  Nx <- nrow(u)
  Ny <- ncol(u)
  
  for (i in 1:Nx) {
    for (j in 1:Ny) {
      # neighbors with Neumann BC (zero-flux)
      uL <- ifelse(i==1, u[i,j], u[i-1,j])
      uR <- ifelse(i==Nx, u[i,j], u[i+1,j])
      uD <- ifelse(j==1, u[i,j], u[i,j-1])
      uU <- ifelse(j==Ny, u[i,j], u[i,j+1])
      
      lap[i,j] <- (uL - 2*u[i,j] + uR)/dx^2 + (uD - 2*u[i,j] + uU)/dy^2
    }
  }
  
  lap
}

# ---- Storage for snapshots ----
snapshots <- list()
snap_times <- c(0, 2, 5, 10, 20)
snap_index <- 1

# ---- Time loop ----
for (n in 1:Nt) {
  
  stim <- alpha*L_ext + beta*x
  
  # reaction term
  Fval <- pgamma(stim, shape=shape, scale=scale)
  reaction <- -x + Fval
  
  # nonlinear diffusion coefficient
  Dmat <- D0 * dgamma(stim, shape=shape, scale=scale)
  
  lap <- laplacian2D(x, dx, dy)
  
  # explicit Euler update
  x_new <- x + dt * (reaction + Dmat * lap)
  
  # keep in [0,1]
  x_new <- pmax(pmin(x_new, 1), 0)
  
  x <- x_new
  
  # save snapshots
  tnow <- n*dt
  if (snap_index <= length(snap_times) && tnow >= snap_times[snap_index]) {
    snapshots[[snap_index]] <- x
    snap_index <- snap_index + 1
  }
}

#save the underlying data
save(snapshots, file = "PDE Model Images/spatial_model_PDE_data.RData")

