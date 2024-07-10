library(data.table)
library(coda)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(stringi)

# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# data tables to save all results
est.data <- data.table()
const.est.data <- data.table()
reassortment.data <- data.table()

#define the number of states
nstates <- 20

# read in the file SIR_simulations.txt that contains the simulated values
simulated <- fread("structuredSIR_simulations.txt", sep="\t")

# get all .log files in /out
logfiles <- list.files("./out", pattern="*.log", full.names=TRUE)

# Process log files in parallel
process_log_file <- function(logfile, simulated, nstates) {
  # Initial processing
  filename <- gsub("./out/", "./master/", logfile)  
  filename <- gsub(".constant.|.ne.|.variable.", ".", filename)  
  lines <- readLines(filename)
  data <- strsplit(lines[2], split="\t")[[1]]
  lins <- str_split(str_replace_all(data[3], "\\[|\\]", ""), ', ')[[1]]
  
  timediff <- as.numeric(strsplit(lins[1], split=":")[[1]][2]) - 
    as.numeric(strsplit(lins[length(lins)], split=":")[[1]][2])
  
  # Set the time points for the parameters
  time_points2 <- seq(0, timediff * 1.1, length.out = 500)
  
  t <- fread(logfile, sep="\t")
  t <- t[round(0.2 * .N):.N]  # remove 20% as burnin
  
  # Calculate ESS for columns 2 to 5
  ess <- sapply(t[, 2:5, with = FALSE], effectiveSize)
  minESS <- min(ess)
  if (minESS < 100) {
    return(NULL)
  }
  
  # Process data based on method type
  method <- strsplit(logfile, "\\.")[[1]][3]
  run <- as.numeric(strsplit(logfile, "\\.|_")[[1]][4])
  
  if (!grepl("constant", method)) {
    prevalence <- sapply(seq_len(length(time_points2) - 1), function(j) {
      if (paste0("reassortment", j) %in% names(t)) {
        list(
          median = median(t[[paste0("reassortment", j)]]),
          lower = quantile(t[[paste0("reassortment", j)]], 0.025),
          upper = quantile(t[[paste0("reassortment", j)]], 0.975),
          time = time_points2[j + 1]
        )
      } else {
        NULL
      }
    })
    prevalence <- rbindlist(prevalence)
    prevalence[, transmission := simulated[run, transmission]]
    prevalence[, c("prevalence", "prevalencel", "prevalenceu") := .(
      median / transmission, lower / transmission, upper / transmission
    )]
    list(method = method, run = run, prevalence = prevalence)
  } else {
    reassortmentRate <- t[["reassortmentRate"]]
    list(
      method = method, run = run,
      prevalence = data.table(
        run = run,
        time = NA,
        prevalence = median(reassortmentRate) / simulated[run, transmission],
        prevalencel = quantile(reassortmentRate, 0.025) / simulated[run, transmission],
        prevalenceu = quantile(reassortmentRate, 0.975) / simulated[run, transmission]
      )
    )
  }
}

# Parallel processing setup
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("simulated", "nstates", "process_log_file"))
clusterEvalQ(cl, {
  library(data.table)
  library(stringi)
  library(coda)
})

results <- parLapply(cl, logfiles, process_log_file, simulated, nstates)
stopCluster(cl)

# Combine results into data tables
for (res in results) {
  if (is.null(res)) next
  if (res$method == "constant") {
    const.est.data <- rbind(const.est.data, res$prevalence)
  } else {
    est.data <- rbind(est.data, res$prevalence)
  }
  reassortment.data <- rbind(reassortment.data, res$reassortment)
}

# Continue with the rest of the script...
