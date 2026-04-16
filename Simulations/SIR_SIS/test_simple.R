library(ggplot2)
library(cowplot)
library(stringr)

# Test just the basic setup
cat("Starting test...\n")

simulated <- read.table("SIS_comp_simulations.txt", header = TRUE, sep = "\t")
cat("Loaded simulations:", nrow(simulated), "\n")

logfiles <- list.files("master", pattern = "*.log", full.names = TRUE)
cat("Found log files:", length(logfiles), "\n")

# Test with just one file
if (length(logfiles) > 0) {
  logfile <- logfiles[1]
  cat("Testing with:", logfile, "\n")

  lines <- readLines(logfile)
  cat("Read", length(lines), "lines\n")

  if (length(lines) >= 2) {
    data <- strsplit(lines[[2]], split = "\t")[[1]]
    cat("Split into", length(data), "parts\n")
  }
}

cat("Test completed successfully!\n")