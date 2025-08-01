# Load necessary packages
library(splines)
library(ggplot2)
library(tidyr)
library(patchwork) # For combining plots
library(dplyr) # For group_by and summarise

# --- Common Parameters ---
p <- 3 # Degree of the B-spline (cubic)
plot_min_x <- 0
plot_max_x <- 10 # Data range from 0 to 10 years
x_vals <- seq(plot_min_x, plot_max_x, length.out = 500) # Dense x-values for smooth curves

# This shift will apply to the simulated data and the fitted spline.
# It ensures visual separation from the basis functions which will stay at y=[0,1].
vertical_shift <- 0.85

# Define colors for basis functions (ensure enough colors for max df = 25)
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#f7b6d2", "#c5b0d5", "#c49c94", "#dbdb8d", "#c7c7c7", "#aec7e8",
            "#ff9896", "#98df8a", "#ffbb78", "#9edae5", "#bcbd22", "#d9d9d9",
            "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d",
            "#666666", "#a6cee3", "#b2df8a") # Total 28+ colors


# --- Simulate Data Points from a Cubic Spline (df=10) with Noise ---
set.seed(42) # For reproducibility
num_data_points <- 500
data_x_simulated <- seq(plot_min_x, plot_max_x, length.out = num_data_points)

# Define parameters for the TRUE underlying cubic spline (df=10)
true_df_for_simulation <- 10
true_degree_for_simulation <- 3 # Cubic spline

# Generate basis matrix for the TRUE spline (evaluated at our data_x_simulated points)
true_basis_matrix_for_data <- bs(data_x_simulated, df = true_df_for_simulation,
                                 degree = true_degree_for_simulation, intercept = TRUE,
                                 Boundary.knots = c(plot_min_x, plot_max_x))

# Define arbitrary "true" coefficients to shape the underlying curve.
# These determine the underlying pattern the spline will try to fit.
# The number of coefficients must match the number of columns in true_basis_matrix_for_data (which is true_df_for_simulation).
# Values chosen to create an interesting shape roughly in the [0, 1] range.
true_coefficients <- c(0.2, 0.8, 0.5, 1.0, 0.3, 0.9, 0.4, 0.7, 0.1, 0.6)

# Calculate the true underlying spline curve
data_y_true <- as.matrix(true_basis_matrix_for_data) %*% true_coefficients

# Add random noise to simulate variability
noise_sd <- 0.08
data_y_noisy <- data_y_true + rnorm(num_data_points, mean = 0, sd = noise_sd)

# Apply the vertical shift to the noisy simulated data for visual separation
sim_data <- data.frame(x = data_x_simulated, y = data_y_noisy + vertical_shift)


# --- Function to Generate and Plot Spline Figure ---
generate_spline_plot <- function(num_basis_functions, plot_title, colors_all, sim_data) {
  
  # Generate B-spline Basis Functions for plotting
  basis_matrix_plot <- bs(x_vals, df = num_basis_functions, degree = p, intercept = TRUE,
                          Boundary.knots = c(plot_min_x, plot_max_x))
  
  # Convert basis matrix for plotting to data frame and reshape
  basis_df_plot <- as.data.frame(basis_matrix_plot)
  colnames(basis_df_plot) <- paste0("B", 1:ncol(basis_df_plot))
  basis_df_plot$x <- x_vals
  basis_long <- pivot_longer(basis_df_plot, cols = starts_with("B"), names_to = "basis_id", values_to = "value")
  
  # Generate basis matrix specifically at the data_x points for regression
  basis_matrix_data <- bs(sim_data$x, df = num_basis_functions, degree = p, intercept = TRUE,
                          Boundary.knots = c(plot_min_x, plot_max_x))
  
  # Perform linear regression to find optimal control points (coefficients)
  fit_model <- lm(sim_data$y ~ basis_matrix_data - 1)
  control_points_optim <- coef(fit_model)
  
  # Calculate the spline curve 's' using the optimized control points
  s_values_fitted <- as.matrix(basis_matrix_plot) %*% control_points_optim
  s_df <- data.frame(x = x_vals, s = s_values_fitted)
  
  # Get the actual internal knots used by bs()
  bs_object_for_knots <- bs(sim_data$x, df = num_basis_functions, degree = p, intercept = TRUE,
                            Boundary.knots = c(plot_min_x, plot_max_x))
  internal_knots_used <- if (!is.null(attr(bs_object_for_knots, "knots"))) as.numeric(attr(bs_object_for_knots, "knots")) else numeric(0)
  # Combine with boundary knots for drawing ticks
  all_knots_to_mark <- unique(c(plot_min_x, internal_knots_used, plot_max_x))
  all_knots_to_mark <- sort(all_knots_to_mark)
  
  
  # --- NEW: Identify and prepare points for red markers on the spline ---
  # Find the x-coordinate where each basis function reaches its peak value
  peak_x_values <- basis_long %>%
    group_by(basis_id) %>%
    summarise(x_peak = x[which.max(value)]) %>%
    ungroup()
  
  # Find the corresponding y-value on the *fitted spline* at these peak x-coordinates
  # This uses the dense s_df to get the spline value at the peak_x_value
  points_on_spline_at_peaks <- data.frame(
    x = peak_x_values$x_peak,
    y = sapply(peak_x_values$x_peak, function(val) {
      s_df$s[which.min(abs(s_df$x - val))] # Find closest x in s_df and get its s value
    })
  )
  # Remove duplicate points if any basis functions peak at the same x (unlikely but possible with discrete x_vals)
  points_on_spline_at_peaks <- unique(points_on_spline_at_peaks)
  
  
  # Determine y-axis limits for the plot to accommodate all elements
  y_axis_max <- max(max(s_df$s), max(sim_data$y), max(basis_long$value)) + 0.1
  y_axis_min <- min(min(s_df$s), min(sim_data$y), min(basis_long$value)) - 0.1
  
  # Adjust y_axis_min to ensure knot ticks are visible below 0
  if (y_axis_min > -0.05) {
    y_axis_min <- -0.05
  }
  
  
  # Plotting
  ggplot() +
    # Simulated Data Points
    geom_point(data = sim_data, aes(x = x, y = y), color = "darkgrey", size = 0.9, alpha = 0.7) +
    
    # Basis functions
    geom_line(data = basis_long, aes(x = x, y = value, color = basis_id), lwd = 1.2) +
    scale_color_manual(values = colors_all) +
    
    # Spline 's'
    geom_line(data = s_df, aes(x = x, y = s), color = "black", lwd = 1.5, linetype = "solid") +
    
    # --- NEW: Add red points on the spline ---
    geom_point(data = points_on_spline_at_peaks, aes(x = x, y = y),
               color = "red", size = 1.1, shape = 19) + # Using solid red circles
    
    
    # Add geom_segment for knot ticks
    geom_segment(aes(x = all_knots_to_mark, xend = all_knots_to_mark,
                     y = -0.025, yend = 0.025), # Tick extending below and above y=0
                 color = "black", lwd = 0.8) +
    
    
    # X-axis setup
    scale_x_continuous(
      breaks = seq(plot_min_x, plot_max_x, by = 1),
      expand = c(0, 0)
    ) +
    # Y-axis setup:
    scale_y_continuous(
      breaks = seq(floor(y_axis_min / 0.2) * 0.2, ceiling(y_axis_max / 0.2) * 0.2, by = 0.2),
      limits = c(y_axis_min, y_axis_max),
      labels = function(x) { ifelse(x < 0, "", x) } # Hide negative labels, effectively starting labels from 0
    ) +
    
    # Theme and labels
    labs(title = plot_title, x = "t", y = "y") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(face="bold"),
      axis.text.y = element_text(size = 10),
      axis.ticks.y = element_line(color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"), # Added Y-axis line
      axis.ticks.x = element_line(color = "black"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

# --- Generate Plots for Different Numbers of Basis Functions ---

# Figure 1: Fewer Basis Functions (df=5)
p1 <- generate_spline_plot(num_basis_functions = 7,
                           plot_title = "a. Small basis",
                           colors_all = colors,
                           sim_data = sim_data)

# Figure 2: More Basis Functions (df=25)
p2 <- generate_spline_plot(num_basis_functions = 19,
                           plot_title = "b. Large basis",
                           colors_all = colors,
                           sim_data = sim_data)

# --- Arrange and Print Figures ---
combined_plot <- p1 | p2 # Stack them horizontally

graphics.off()
grDevices::cairo_pdf("res/bs_ex.pdf", height = 4.5, width = 10)
combined_plot
dev.off()
