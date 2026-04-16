library(ggplot2)
library(cowplot)
library(stringr)

# Clear workspace
rm(list = ls())

# Set the directory to the directory of the file
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)
# Use current working directory instead

# Read simulation parameters
simulated <- read.table("SIS_comp_simulations.txt", header = TRUE, sep = "\t")

# Get all .log files in master/
logfiles <- list.files("master", pattern = "*.log", full.names = TRUE)

# Helper: parse a single log file and compute metrics
parse_log <- function(logfile, population_size) {
  lines <- readLines(logfile)
  data <- strsplit(lines[[2]], split = "\t")[[1]]

  sir <- str_split(str_replace_all(data[[2]], "\\[|\\]", ""), ", ")[[1]]

  I <- 1
  time_I <- 0
  co_inf_count <- 0
  co_inf_times <- numeric()

  for (j in seq_along(sir)) {
    tmp <- str_split(sir[j], ":")[[1]]
    event <- tmp[1]
    t_event <- as.numeric(tmp[2])
    if (event == "0") {
      I <- c(I, I[length(I)])
      time_I <- c(time_I, t_event)
      co_inf_count <- co_inf_count + 1
      co_inf_times <- c(co_inf_times, t_event)
    } else if (event == "1") {
      I <- c(I, I[length(I)] + 1)
      time_I <- c(time_I, t_event)
    } else if (event == "2") {
      I <- c(I, I[length(I)] - 1)
      time_I <- c(time_I, t_event)
    } else if (event == "3") {
      I <- c(I, I[length(I)])
      time_I <- c(time_I, t_event)
    }
  }

  # Peak prevalence and its time
  peak_idx <- which.max(I)
  peak_prev <- max(I) / population_size
  peak_prev_time <- time_I[peak_idx]

  # Peak co-infection time via kernel density
  if (length(co_inf_times) >= 3) {
    d <- density(co_inf_times)
    peak_coinf_time <- d$x[which.max(d$y)]
  } else {
    peak_coinf_time <- NA
  }

  # Time-weighted average prevalence
  dt <- diff(time_I)
  if (length(dt) > 0) {
    avg_prev <- sum(I[-length(I)] * dt) / sum(dt)
  } else {
    avg_prev <- I[1]
  }

  # Normalized co-infections
  co_inf_per_prev <- if (avg_prev > 0) co_inf_count / avg_prev else NA

  list(
    peak_prevalence = peak_prev,
    peak_prevalence_time = peak_prev_time,
    peak_coinf_time = peak_coinf_time,
    total_co_infections = co_inf_count,
    average_prevalence = avg_prev,
    co_inf_per_prevalence = co_inf_per_prev,
    I = I,
    time_I = time_I
  )
}

# Parse all log files and build results data frame
results <- data.frame()
trajectories <- data.frame()

for (logfile in logfiles) {
  # Extract model type and run number from filename
  basename_f <- basename(logfile)
  if (grepl("^SIR_simulations_", basename_f)) {
    model <- "SIR"
  } else if (grepl("^SIS_simulations_", basename_f)) {
    model <- "SIS"
  } else {
    next
  }
  run <- as.numeric(str_extract(basename_f, "\\d+"))

  pop_size <- simulated[run, "population_size"]
  transmission <- simulated[run, "transmission"]
  R0 <- transmission  # since recovery = 1

  parsed <- tryCatch(parse_log(logfile, pop_size), error = function(e) NULL)
  if (is.null(parsed)) next

  results <- rbind(results, data.frame(
    run = run,
    model = model,
    R0 = R0,
    population_size = pop_size,
    peak_prevalence = parsed$peak_prevalence,
    peak_prevalence_time = parsed$peak_prevalence_time,
    peak_coinf_time = parsed$peak_coinf_time,
    total_co_infections = parsed$total_co_infections,
    average_prevalence = parsed$average_prevalence,
    co_inf_per_prevalence = parsed$co_inf_per_prevalence
  ))

  # Store downsampled trajectory
  n <- length(parsed$time_I)
  idx <- unique(c(1, seq(1, n, length.out = min(200, n)), n))
  trajectories <- rbind(trajectories, data.frame(
    run = run,
    model = model,
    R0 = R0,
    time = parsed$time_I[idx],
    I = parsed$I[idx] / pop_size
  ))
}

# ===== Figure 1: Raw comparison (3 panels) =====

# Panel A: Prevalence trajectories for 6 representative R0 values
R0_vals <- sort(unique(results$R0))
target_R0s <- quantile(R0_vals, probs = c(0.05, 0.25, 0.45, 0.55, 0.75, 0.95))
# Pick the closest actual R0 for each target
rep_runs <- sapply(target_R0s, function(r) {
  which.min(abs(results$R0[results$model == "SIR"] - r))
})
rep_run_ids <- unique(results$run[results$model == "SIR"])[rep_runs]

traj_sub <- trajectories[trajectories$run %in% rep_run_ids, ]
# Create R0 label for faceting
traj_sub$R0_label <- sprintf("R0 = %.1f", traj_sub$R0)
# Normalize time: shift so peak is at t=0 for each run/model
for (r in unique(traj_sub$run)) {
  for (m in c("SIR", "SIS")) {
    mask <- traj_sub$run == r & traj_sub$model == m
    if (sum(mask) > 0) {
      peak_time <- traj_sub$time[mask][which.max(traj_sub$I[mask])]
      traj_sub$time[mask] <- traj_sub$time[mask] - peak_time
    }
  }
}

pA <- ggplot(traj_sub, aes(x = time, y = I, color = model)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~R0_label, scales = "free", nrow = 1) +
  scale_color_manual(values = c("SIR" = "#0072B2", "SIS" = "#D55E00")) +
  labs(x = "Time (centered on peak)", y = "Prevalence (I/N)", color = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Panel B: Peak prevalence vs R0
pB <- ggplot(results, aes(x = R0, y = peak_prevalence, color = model)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("SIR" = "#0072B2", "SIS" = "#D55E00")) +
  labs(x = "R0", y = "Peak prevalence (I/N)", color = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Panel C: Co-infections / avg prevalence vs R0
pC <- ggplot(results, aes(x = R0, y = co_inf_per_prevalence, color = model)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = c("SIR" = "#0072B2", "SIS" = "#D55E00")) +
  labs(x = "R0", y = "Co-infections / avg prevalence", color = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom")

fig1 <- plot_grid(
  pA,
  plot_grid(pB, pC, ncol = 2, labels = c("B", "C")),
  nrow = 2, labels = c("A", ""), rel_heights = c(1, 1)
)

ggsave("../../Figures/SIR_SIS_comparison.pdf", fig1, width = 12, height = 8)

# ===== Figure 1b: Reassortment rate vs prevalence over time =====

# Calculate reassortment rate over time for representative runs
reassort_data <- data.frame()

for (logfile in logfiles) {
  basename_f <- basename(logfile)
  if (grepl("^SIR_simulations_", basename_f)) {
    model <- "SIR"
  } else if (grepl("^SIS_simulations_", basename_f)) {
    model <- "SIS"
  } else {
    next
  }
  run <- as.numeric(str_extract(basename_f, "\\d+"))

  # Only process representative runs used in panel A
  if (!run %in% rep_run_ids) next

  pop_size <- simulated[run, "population_size"]

  # Parse log file to get detailed time series
  lines <- readLines(logfile)
  data <- strsplit(lines[[2]], split = "\t")[[1]]
  sir <- str_split(str_replace_all(data[[2]], "\\[|\\]", ""), ", ")[[1]]

  I <- 1
  time_I <- 0
  co_inf_times <- numeric()

  for (j in seq_along(sir)) {
    tmp <- str_split(sir[j], ":")[[1]]
    event <- tmp[1]
    t_event <- as.numeric(tmp[2])
    if (event == "0") {
      I <- c(I, I[length(I)])
      time_I <- c(time_I, t_event)
      co_inf_times <- c(co_inf_times, t_event)
    } else if (event == "1") {
      I <- c(I, I[length(I)] + 1)
      time_I <- c(time_I, t_event)
    } else if (event == "2") {
      I <- c(I, I[length(I)] - 1)
      time_I <- c(time_I, t_event)
    } else if (event == "3") {
      I <- c(I, I[length(I)])
      time_I <- c(time_I, t_event)
    }
  }

  # Calculate reassortment rate using 200 equally spaced time points
  if (length(co_inf_times) >= 3) {
    time_points <- seq(0, max(time_I), length.out = 200)  # 200 equally spaced points
    reassort_rates <- numeric(length(time_points) - 1)
    window_times <- numeric(length(time_points) - 1)
    window_prevalence <- numeric(length(time_points) - 1)

    for (i in 1:(length(time_points) - 1)) {
      t_start <- time_points[i]
      t_end <- time_points[i + 1]
      window_times[i] <- (t_start + t_end) / 2

      # Count co-infections in this time window
      coinf_in_window <- sum(co_inf_times >= t_start & co_inf_times < t_end)

      # Calculate average number of infected individuals (not prevalence) in this time window
      time_mask <- time_I >= t_start & time_I <= t_end
      mean_infected <- 0  # Initialize to avoid scoping issues

      if (sum(time_mask) > 0) {
        mean_infected <- mean(I[time_mask])
        window_prevalence[i] <- mean_infected / pop_size
      } else {
        # Interpolate number of infected individuals
        if (t_start <= max(time_I)) {
          prev_idx <- max(which(time_I <= t_start))
          next_idx <- min(which(time_I >= t_end), length(time_I))
          mean_infected <- mean(I[c(prev_idx, next_idx)])
          window_prevalence[i] <- mean_infected / pop_size
        } else {
          mean_infected <- tail(I, 1)
          window_prevalence[i] <- mean_infected / pop_size
        }
      }

      # Calculate per-individual reassortment rate
      if (mean_infected > 0) {
        reassort_rates[i] <- (coinf_in_window / (t_end - t_start)) / mean_infected
      } else {
        reassort_rates[i] <- 0  # No reassortment possible when no infected individuals
      }
    }

    # Normalize time relative to peak prevalence
    peak_time <- time_I[which.max(I)]

    reassort_data <- rbind(reassort_data, data.frame(
      run = run,
      model = model,
      R0 = simulated[run, "transmission"],
      time = window_times - peak_time,
      reassortment_rate = reassort_rates,
      prevalence = window_prevalence,
      R0_label = sprintf("R0 = %.1f", simulated[run, "transmission"])
    ))
  }
}

# Create combined reassortment plot with SIR prevalence and both SIR/SIS reassortment rates
library(scales)

# Prepare data for plotting
sir_data <- reassort_data[reassort_data$model == "SIR", ]
sis_data <- reassort_data[reassort_data$model == "SIS", ]

if (nrow(sir_data) > 0 && nrow(sis_data) > 0) {

  # Create R0 labels for subplot titles (round to 1 decimal)
  sir_data$R0_rounded <- round(sir_data$R0, 1)
  sir_data$R0_title <- paste("R0 =", sir_data$R0_rounded)

  sis_data$R0_rounded <- round(sis_data$R0, 1)
  sis_data$R0_title <- paste("R0 =", sis_data$R0_rounded)

  # Normalize each metric to its OWN maximum within each R0 group
  sir_data <- do.call(rbind, lapply(split(sir_data, sir_data$R0_title), function(df) {
    df$prevalence_norm <- df$prevalence / max(df$prevalence, na.rm = TRUE)
    df$reassortment_rate_norm <- df$reassortment_rate / max(df$reassortment_rate, na.rm = TRUE)
    return(df)
  }))

  sis_data <- do.call(rbind, lapply(split(sis_data, sis_data$R0_title), function(df) {
    df$reassortment_rate_norm <- df$reassortment_rate / max(df$reassortment_rate, na.rm = TRUE)
    return(df)
  }))

  # Create the combined plot using separate geom_line calls
  p_combined <- ggplot() +
    # SIR prevalence
    geom_line(data = sir_data, aes(x = time, y = prevalence_norm, color = "SIR Prevalence"),
              linewidth = 1.2, alpha = 0.9) +
    # SIR co-infection rate (per individual)
    geom_line(data = sir_data, aes(x = time, y = reassortment_rate_norm, color = "SIR Co-infection rate (per individual)"),
              linewidth = 1.2, alpha = 0.9) +
    # SIS co-infection rate (per individual)
    geom_line(data = sis_data, aes(x = time, y = reassortment_rate_norm, color = "SIS Co-infection rate (per individual)"),
              linewidth = 1.2, alpha = 0.9) +
    facet_wrap(~ factor(R0_title, levels = unique(sir_data$R0_title[order(sir_data$R0_rounded)])),
               scales = "free_x", nrow = 2) +
    scale_color_manual(
      values = c(
        "SIR Prevalence" = "#0072B2",           # Blue (matching SIR style)
        "SIR Co-infection rate (per individual)" = "#D55E00",    # Orange (matching SIR style)
        "SIS Co-infection rate (per individual)" = "#009E73"     # Green
      ),
      name = ""
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.0", "0.25", "0.50", "0.75", "1.0")
    ) +
    labs(
      x = "Time (centered on epidemic peak)",
      y = "Normalized intensity"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      strip.text = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

  # Save the combined plot
  ggsave("../../Figures/SIR_SIS_reassortment_comparison.pdf", p_combined, width = 12, height = 8)
}

# ===== Figure 2: Ratio comparison vs R0 =====

# Merge paired data
sir_data <- results[results$model == "SIR", ]
sis_data <- results[results$model == "SIS", ]
paired <- merge(sir_data, sis_data, by = "run", suffixes = c("_SIR", "_SIS"))

paired$peak_ratio <- paired$peak_prevalence_SIS / paired$peak_prevalence_SIR
paired$coinf_ratio <- paired$co_inf_per_prevalence_SIS / paired$co_inf_per_prevalence_SIR
# Express time differences in generations (multiply by recovery rate)
recovery_rate <- simulated[paired$run, "recovery"]
paired$peak_prev_time_diff <- (paired$peak_prevalence_time_SIS - paired$peak_prevalence_time_SIR) * recovery_rate
paired$peak_coinf_time_diff <- (paired$peak_coinf_time_SIS - paired$peak_coinf_time_SIR) * recovery_rate

# Panel A: Normalized co-infection ratio
p2A <- ggplot(paired, aes(x = R0_SIR, y = coinf_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 2, color = "#CC79A7") +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, color = "#CC79A7", fill = "#CC79A7", alpha = 0.2) +
  labs(x = "R0", y = "SIS / SIR ratio") +
  theme_minimal()

# Panel B: Difference in peak co-infection time (SIS - SIR)
p2B <- ggplot(paired[!is.na(paired$peak_coinf_time_diff), ],
              aes(x = R0_SIR, y = peak_coinf_time_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 2, color = "#D55E00") +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, color = "#D55E00", fill = "#D55E00", alpha = 0.2) +
  labs(x = "R0", y = "Time difference [generations]") +
  theme_minimal()

fig2 <- plot_grid(p2A, p2B, ncol = 2, labels = c("A", "B"))

ggsave("../../Figures/SIR_SIS_ratios.pdf", fig2, width = 10, height = 4)
