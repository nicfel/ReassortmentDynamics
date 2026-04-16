library(stringr)
library(coda)
library(ggplot2)
# add library for color palette
library(RColorBrewer)


# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# date frame to save all results
est.data = data.frame()
# data frame to save the reassortment events
reassortment.data = data.frame()



# read in the file SIR_simulations.txt that contains the simulated values
simulated = read.table("SIR_simulations.txt", header=TRUE, sep="\t")

# get the true values for the prevalence by readin in the files in /master
true_data = data.frame()
co_inf = data.frame()
inf = data.frame()
# get all .log files in .master
logfiles <- list.files("./out", pattern="*.log", full.names=TRUE)

# loop over all files
for (i in seq(1, length(logfiles))){
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][4])

  # read in the file line by line
  lines = readLines(logfiles[i])
  for (k in seq(2, 2)){#length(lines))){
    data = strsplit(lines[[k]], split="\t")[[1]]

    # Initialize lists to hold the results
    I_list <- list()
    time_I_list <- list()
    co_inf_times_list <- list()
  
    # Replace [ and ] then split the SIR events into a vector of strings
    sir <- str_split(str_replace_all(data[[2]], "\\[|\\]", ""), ', ')[[1]]
    # Do the same for lineages
    lins <- str_split(str_replace_all(data[[3]], "\\[|\\]", ""), ', ')[[1]]
    
    # get the most recent sampled individual
    mrsi = strsplit(lins[[1]], split=":")[[1]][[2]]
    # add the mrsi to the est.data
    # est.data[est.data$run==run, "mrsi"] = as.numeric(mrsi)
    # 
    # Compute the number of infected individuals over time
    I <- 1
    time_I <- 0
    # Keep also track of the co-infection event times
    co_inf_times <- numeric()
    inf_times <- numeric()
  
    for (j in seq_along(sir)) {
        tmp <- str_split(sir[j], ':')[[1]]
        if (tmp[1] == '0') {
            I <- c(I, I[length(I)])
            time_I <- c(time_I, as.numeric(tmp[2]))
            co_inf_times <- c(co_inf_times, time_I[length(time_I)])
        } else if (tmp[1] == '1') {
            I <- c(I, I[length(I)] + 1)
            time_I <- c(time_I, as.numeric(tmp[2]))
            inf_times <- c(inf_times, time_I[length(time_I)])
        } else if (tmp[1] == '2') {
            I <- c(I, I[length(I)] - 1)
            time_I <- c(time_I, as.numeric(tmp[2]))
        } else if (tmp[1] == '3') {
            I <- c(I, I[length(I)])
            time_I <- c(time_I, as.numeric(tmp[2]))
        }
    }
  
    # get the indec for which the I is max
    max_index <- which.max(I)
    co_inf_times <- co_inf_times - time_I[max_index]
    inf_times <- inf_times - time_I[max_index]
    time_I <- time_I-time_I[max_index]

    # get the population size
    popSize = simulated[run, "population_size"]
    
    # take 100 points from the time_I vector
    index = seq(1, length(time_I), length.out=500)
    
    if (grepl("SIS", logfiles[i])){
      # if the file is a SIS file, then we need to take the first 100 points
      method = "transmission rate limited"
    }else{
      method = "susceptible depletion"
    }
    
    # add the data to the data frame
    true_data = rbind(true_data, data.frame(run=k, time=time_I[index], prevalence=I[index], method))
    co_inf = rbind(co_inf, data.frame(run=k, time=co_inf_times, method))
    inf = rbind(inf, data.frame(run=k, time=inf_times, method))
  }
}

# Estimate density of co-infection events over time for each method
# Create time grid for density estimation
time_range <- range(co_inf$time)
time_grid <- seq(time_range[1], time_range[2], length.out = 100)

# Calculate density estimates for each method
co_inf_density <- data.frame()
for (method_type in unique(co_inf$method)) {
  co_inf_subset <- co_inf[co_inf$method == method_type, ]
  
  if (nrow(co_inf_subset) > 0) {
    # Estimate density of co-infection events
    dens <- density(co_inf_subset$time, from = time_range[1], to = time_range[2], n = 100)
    co_inf_density <- rbind(co_inf_density, 
                            data.frame(time = dens$x, 
                                       density = dens$y,
                                       method = method_type))
  }
}

# Create prevalence density for comparison (average prevalence over time)
prevalence_avg <- aggregate(prevalence ~ time + method, data = true_data, FUN = mean)
prevalence_avg <- aggregate(prevalence ~ time + method, data = true_data, FUN = mean)

# Create the plots without individual legends
p1 = ggplot(true_data, aes(x=time, y=prevalence, color=method))+
  geom_density(data=co_inf, aes(y=..count.., color=method), fill=NA) +
  geom_density(data=inf, aes(y=..count.., color=method), fill=NA) +
  theme_minimal() +
  coord_cartesian(xlim= c(-3, 5)) +
  xlab("Time relative to peak") +
  ylab("Number of Co-Infection Events") +
  theme(legend.position = "none")  # Remove legend
plot(p1)

p2 = ggplot(true_data, aes(x=time, y=prevalence, color=method, group=interaction(method, run)))+
  geom_line() +
  theme_minimal() +
  coord_cartesian(xlim= c(-3, 5)) +
  xlab("Time relative to peak") +
  ylab("Prevalence") +
  theme(legend.position = "none")  # Remove legend

time_range <- range(co_inf$time)
time_breaks <- seq(time_range[1], time_range[2], length.out = 51)  # 100 bins
time_centers <- (time_breaks[-1] + time_breaks[-length(time_breaks)]) / 2

co_inf_normalized <- data.frame()

for (method_type in unique(co_inf$method)) {
  co_inf_subset <- co_inf[co_inf$method == method_type, ]
  true_data_subset <- true_data[true_data$method == method_type, ]
  
  if (nrow(co_inf_subset) > 0 & nrow(true_data_subset) > 0) {
    # Count co-infection events in each time bin
    co_inf_counts <- hist(co_inf_subset$time, breaks = time_breaks, plot = FALSE)$counts
    
    # Calculate average prevalence at the same time points
    avg_prevalence <- aggregate(prevalence ~ time, data = true_data_subset, FUN = mean)
    prevalence_interp <- approx(avg_prevalence$time, avg_prevalence$prevalence, 
                                xout = time_centers, rule = 2)$y
    
    # Normalize counts by prevalence (add small constant to avoid division by zero)
    normalized_counts <- co_inf_counts / (prevalence_interp + 1e-6)
    
    co_inf_normalized <- rbind(co_inf_normalized, 
                               data.frame(time = time_centers, 
                                          normalized_counts = normalized_counts,
                                          method = method_type))
  }
}

# save co_inf_normalized to a file
write.table(co_inf_normalized, file = "co_inf_normalized.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create the normalized counts plot for p3
p3 = ggplot() +
  geom_line(data = co_inf_normalized, aes(x = time, y = normalized_counts, 
                                          color = method), size = 1) +
  theme_minimal() +
  xlab("Time relative to peak") +
  ylab("Normalized Co-infection Counts") +
  coord_cartesian(xlim = c(-3, 5)) +
  labs(title = "Co-infection counts normalized by prevalence") +
  guides(color = guide_legend(title = "Method")) +
  theme(legend.position = "none")

# Extract legend from one of the plots
legend <- get_legend(
  p3 + 
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Method"))
)

# Combine plots without legends
plots_combined <- plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 1)

# Add the shared legend at the bottom
combined_plot <- plot_grid(plots_combined, legend, ncol = 1, rel_heights = c(1, 0.1))

plot(combined_plot)
ggsave("combined_plot.png", plot = combined_plot, width = 10, height = 15)

# Print summary of co-infection timing
cat("Summary of Co-infection Event Timing by Method:\n")
summary_timing <- aggregate(time ~ method, data = co_inf, 
                            FUN = function(x) c(mean = mean(x), median = median(x), 
                                                sd = sd(x), min = min(x), max = max(x)))
print(summary_timing)

# Save the final plot (fixed the variable name)
ggsave(filename = "../../Figures/SIR_SIR.pdf", plot = combined_plot,
       width = 12, height = 8)