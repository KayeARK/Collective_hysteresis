# Morning Chorus Model: Catastrophe Theory Analysis of Bird Singing Behavior

This repository contains R code for analyzing collective bird singing behavior using mathematical models based on catastrophe theory and cusp manifold analysis. The project examines how environmental light levels and social coupling influence the onset of collective singing in bird populations.

## R Environment Requirements

**R Version:** R ≥ 4.0.0 (tested with R 4.3.x)

**Required Packages:**
```r
install.packages(c(
  "ggplot2",      # Core plotting and visualization
  "dplyr",        # Data manipulation and filtering
  "plotly",       # Interactive 3D plots
  "tidyr",        # Data reshaping
  "RColorBrewer", # Color palettes for scientific plots
  "reshape2",     # Data transformation
  "viridis",      # Perceptually uniform color scales
  "grid",         # Grid graphics system
  "htmlwidgets"   # Save interactive plots
))
```

**Optional Dependencies:**
- **pdflatex**: For compiling figure outputs (if generating publication figures)
- **Pandoc**: For generating documentation from R Markdown (if applicable)

## Computational Requirements

### Runtime Estimates
- **Basic bifurcation analysis** (`number_of_fps.R`): ~2-5 seconds
- **Parameter variation studies** (`*_vary_params_*.R`): ~10-30 seconds each
- **3D manifold visualization** (`catstrophe_manifold_plot.R`): ~2-5 minutes
- **Agent-based simulations** (`Model_gamma_individual_*.r`): ~5-15 minutes each
- **Spatial models** (`Spatial_model_*.R`): ~10-30 minutes each

**Total runtime for all analyses:** ~1-2 hours

### System Requirements
- **RAM:** Minimum 4GB, recommended 8GB+ (for large agent-based simulations)
- **CPU:** Multi-core processor recommended for spatial models
- **Storage:** ~100MB for code and generated figures

## Project Structure

```
Morning chorus/
 Code/
   ├── Figure 1/           # Basic gamma distribution plot for schematic
   │   └── gamma_distribution_plot.R
   ├── Figure 2/           # Core bifurcation analysis & agent models
   │   ├── number_of_fps.R                    # Main fold curve analysis
   │   ├── catstrophe_manifold_plot.R         # 3D cusp manifold
   │   └── Model_gamma_individual_b[3-5].r    # Agent-based simulations
   ├── Figure 3/           # Parameter sensitivity (gamma parameters)
   │   ├── number_of_fps_vary_params_mean.R   # Vary distribution mean
   │   ├── number_of_fps_vary_params_scale.R  # Vary scale parameter
   │   └── number_of_fps_vary_params_shape.R  # Vary shape parameter
   ├── Figure 4/           # Alternative distributions
   │   ├── number_of_fps_vary_params_beta.R     # Beta distributions
   │   ├── number_of_fps_vary_params_lognormal.R # Log-normal distributions
   │   └── number_of_fps_vary_params_weibull.R  # Weibull distributions
   └── Figure 5/           # Spatial models
       ├── Spatial_model_individual.R         # Individual-based spatial model
       ├── Spatial_model_PDE.R                # PDE-based spatial model
       └── Figure generation.R                # Combine spatial outputs

```

## Quick Start Guide

### 1. Basic Analysis
Start with the core bifurcation analysis:
```r
source("Code/Figure 2/number_of_fps.R")
```
This generates the fundamental fold curve showing regions of 1 vs 3 fixed points.

### 2. 3D Visualization
Create interactive cusp manifold:
```r
source("Code/Figure 2/catstrophe_manifold_plot.R")
```
Outputs `cusp_catastrophe_gamma.html` for interactive exploration.

### 3. Parameter Studies
Explore sensitivity to gamma distribution parameters:
```r
source("Code/Figure 3/number_of_fps_vary_params_shape.R")
source("Code/Figure 3/number_of_fps_vary_params_scale.R")  
source("Code/Figure 3/number_of_fps_vary_params_mean.R")
```

### 4. Agent-Based Validation
Run population simulations (computationally intensive):
```r
source("Code/Figure 2/Model_gamma_individual_b4.r")  # β = 4
```

### 5. Spatial Dynamics
Generate spatial propagation patterns:
```r
source("Code/Figure 5/Spatial_model_individual.R")
```

## Core Mathematical Framework

### Key Variables
- **L**: Light level (environmental drive, range 0-1.5)
- **β**: Social coupling strength (interaction parameter, range 0-6)
- **α**: Light sensitivity coefficient (typically 1.0)
- **shape, scale**: Gamma distribution parameters for individual thresholds

### Central Functions
All analyses use consistent mathematical transformations:
```r
f_gamma(z) = dgamma(z, shape, scale)              # Probability density
F_gamma(z) = pgamma(z, shape, scale)              # Cumulative distribution
beta_of_z(z) = 1 / f_gamma(z)                    # Social coupling curve
L_of_z(z) = (z - F_gamma(z)/f_gamma(z)) / alpha  # Light level curve
```

## Output Files

### Generated Figures
- **PDF outputs**: High-resolution figures for publication
- **Interactive HTML**: 3D manifold visualizations
- **Data files**: `.RData` files containing simulation results for reanalysis

### File Naming Convention
- `analytic_*`: Theoretical/mathematical analysis
- `*_vary_params_*`: Parameter sensitivity studies  
- `spatial_*`: Spatial model outputs
- `hysteresis_*`: Agent-based simulation results

## Troubleshooting

### Common Issues
1. **Memory errors in spatial models**: Reduce grid resolution (`Nx`, `Ny` parameters)
2. **Slow agent simulations**: Decrease population size (`N` parameter) or time steps (`T_max`)
3. **Missing interactive plots**: Ensure `plotly` and `htmlwidgets` packages installed
4. **PDF generation errors**: Check that output directories exist

### Performance Optimization
- Run computationally expensive scripts (`Model_*`, `Spatial_*`) separately
- Pre-computed results stored in `.RData` files can be loaded to skip lengthy calculations
- Adjust resolution parameters (`length.out`, `N`, grid sizes) based on available computational resources

## Citation

If using this code for research, please cite the associated manuscript:

*[Citation details to be added upon publication]*

## License

This project is licensed under the MIT License - see `LICENSE.txt` for details.