# Required packages
library(plotly)
library(dplyr)

# Gamma CDF parameters (choose shape > 1 for inflection)
shape <- 4
scale <- 0.8

F <- function(z) pgamma(z, shape = shape, scale = scale)

alpha <- 1   # light sensitivity

L_vals    <- seq(0, 1.5, length.out = 400)     # light level
beta_vals <- seq(0, 6, length.out = 400)     # social coupling
#beta_vals <- 4
x_grid    <- seq(0, 1, length.out = 400)    # possible steady states

find_roots <- function(L, beta) {
  g <- function(x) -x + F(alpha * L + beta * x)
  vals <- g(x_grid)
  idx <- which(diff(sign(vals)) != 0)

  roots <- sapply(idx, function(i) {
    uniroot(g, c(x_grid[i], x_grid[i+1]))$root
  })

  roots
}

stability <- function(x, L, beta) {
  dF <- dgamma(alpha * L + beta * x,
               shape = shape, scale = scale)
  deriv <- -1 + beta * dF
  ifelse(deriv < 0, "stable", "unstable")
}

data <- list()

for (L in L_vals) {
  for (beta in beta_vals) {
    roots <- find_roots(L, beta)
    if (length(roots) > 0) {
      for (x in roots) {
        data[[length(data) + 1]] <- data.frame(
          L = L,
          beta = beta,
          x = x,
          stability = stability(x, L, beta)
        )
      }
    }
  }
}

df <- bind_rows(data)

# Find the actual beta values closest to 3, 4, and 5
target_betas <- c(3, 4, 5)
actual_beta_values <- c()

for (target in target_betas) {
  closest_beta <- beta_vals[which.min(abs(beta_vals - target))]
  actual_beta_values <- c(actual_beta_values, closest_beta)
}

cat("Target betas:", target_betas, "\n")
cat("Actual betas:", actual_beta_values, "\n")

# Add a column to identify special beta values
df$color_group <- ifelse(
  df$beta %in% actual_beta_values,
  "special_beta",
  ifelse(df$stability == "stable", "stable", "unstable")
)

# Debug: Check how many points we have in each category
cat("Points by color group:\n")
print(table(df$color_group))
cat("Sample of special beta points:\n")
special_points <- df[df$color_group == "special_beta", ]
if(nrow(special_points) > 0) {
  print(head(special_points[c("L", "beta", "x", "stability")], 10))
} else {
  cat("No special beta points found!\n")
}

p <- plot_ly()

p <- p %>%
  add_markers(
    data = df %>% filter(color_group == "stable"),
    x = ~L, y = ~beta, z = ~x,
    marker = list(size = 2, color = "#1f77b4"),
    name = "Stable"
  ) %>%
  add_markers(
    data = df %>% filter(color_group == "unstable"),
    x = ~L, y = ~beta, z = ~x,
    marker = list(size = 2, color = "#d62728"),
    name = "Unstable"
  ) %>%
  add_markers(
    data = df %>% filter(color_group == "special_beta"),
    x = ~L, y = ~beta, z = ~x,
    marker = list(size = 3, color = "black", symbol = "circle"),
    name = "Beta = 3, 4, 5"
  )

p <- p %>%
  layout(
    scene = list(
      xaxis = list(title = "Light level, L", showgrid = FALSE),
      yaxis = list(title = "Social coupling, β", showgrid = FALSE),
      zaxis = list(title = "Steady state, x", showgrid = FALSE)
    )
  )

p

# To save the plot as an HTML file
htmlwidgets::saveWidget(p, "cusp_catastrophe_gamma.html")