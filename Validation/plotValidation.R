library(stringr)
library(coda)
library(ggplot2)
library(RColorBrewer)
library(colorblindr)
library(patchwork)

# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Define simulation names
sim_names <- c("constant", "skygrowth_1", "skygrowth_2")

# Initialize lists to store data
sim_data <- list()
inf_data <- list()

# Read simulation log files (simulated under the prior)
for (sim_name in sim_names) {
  # Handle naming differences: constant vs skygrowth files
  if (sim_name == "constant") {
    sim_file <- paste0("simulate_", sim_name, ".log")
  } else {
    sim_file <- paste0("simulate_all_", sim_name, ".log")
  }
  if (file.exists(sim_file)) {
    sim_data[[sim_name]] <- read.table(sim_file, header=TRUE, sep="\t")
    cat("Read", sim_file, "-", nrow(sim_data[[sim_name]]), "samples\n")
  } else {
    cat("Warning: File", sim_file, "not found\n")
  }
}

# Read inference log files (sampled under the prior)
for (sim_name in sim_names) {
  inf_file <- paste0("infer_all_", sim_name, ".log")
  if (file.exists(inf_file)) {
    inf_data[[sim_name]] <- read.table(inf_file, header=TRUE, sep="\t")
    cat("Read", inf_file, "-", nrow(inf_data[[sim_name]]), "samples\n")
  } else {
    cat("Warning: File", inf_file, "not found\n")
  }
}

# Prepare data for plotting
plot_data <- data.frame()

for (sim_name in sim_names) {
  if (sim_name %in% names(sim_data) && sim_name %in% names(inf_data)) {
    # Network length
    sim_length <- sim_data[[sim_name]]$network.totalLength
    inf_length <- inf_data[[sim_name]]$network.totalLength
    
    # Reassortment node count
    sim_reassort <- sim_data[[sim_name]]$network.reassortmentNodeCount
    inf_reassort <- inf_data[[sim_name]]$network.reassortmentNodeCount
    
    # Add simulation data
    plot_data <- rbind(plot_data,
      data.frame(
        Simulation = sim_name,
        Metric = "Network Length",
        Value = sim_length,
        Type = "Simulated"
      ),
      data.frame(
        Simulation = sim_name,
        Metric = "Network Length",
        Value = inf_length,
        Type = "Sampled"
      ),
      data.frame(
        Simulation = sim_name,
        Metric = "Reassortment Node Count",
        Value = sim_reassort,
        Type = "Simulated"
      ),
      data.frame(
        Simulation = sim_name,
        Metric = "Reassortment Node Count",
        Value = inf_reassort,
        Type = "Sampled"
      )
    )
  }
}

# Clean up simulation names for better labels
plot_data$Simulation <- factor(plot_data$Simulation,
  levels = c("constant", "skygrowth_1", "skygrowth_2"),
  labels = c("Constant", "Skygrowth 1", "Skygrowth 2")
)

plot_data$Type <- factor(plot_data$Type, levels = c("Simulated", "Sampled"))

# Separate data for network length and reassortment node count
length_data <- plot_data[plot_data$Metric == "Network Length", ]
reassort_data <- plot_data[plot_data$Metric == "Reassortment Node Count", ]

# Get range for reassortment node count to set x-axis limits
reassort_max <- max(reassort_data$Value, na.rm = TRUE)
reassort_min <- min(reassort_data$Value, na.rm = TRUE)

# Create plot for network length (density)
p_length <- ggplot(length_data, aes(x = Value, fill = Type, color = Type)) +
  geom_density(alpha = 0.3, adjust = 1.5) +
  facet_grid(. ~ Simulation, scales = "free") +
  scale_fill_manual(values = c("Simulated" = "#E69F00", "Sampled" = "#56B4E9")) +
  scale_color_manual(values = c("Simulated" = "#E69F00", "Sampled" = "#56B4E9")) +
  labs(
    x = "Network Length",
    y = "Density",
    fill = "Type",
    color = "Type"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.position = "none"
  )

# Create plot for reassortment node count (histogram)
p_reassort <- ggplot(reassort_data, aes(x = Value, fill = Type, color = Type)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  facet_grid(. ~ Simulation, scales = "free") +
  scale_fill_manual(values = c("Simulated" = "#E69F00", "Sampled" = "#56B4E9")) +
  scale_color_manual(values = c("Simulated" = "#E69F00", "Sampled" = "#56B4E9")) +
  scale_x_continuous(limits = c(reassort_min, reassort_max)) +
  labs(
    x = "Reassortment Node Count",
    y = "Density",
    fill = "Type",
    color = "Type",
    title = ""
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# Combine plots
p <- p_length / p_reassort + plot_layout(heights = c(1, 1.1))

# Save the plot
ggsave("../Figures/Validation.pdf", p, width = 10, height = 6)

cat("\nPlot saved to ../Figures/Validation.pdf\n")
