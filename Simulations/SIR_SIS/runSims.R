# Clear workspace
rm(list = ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Make a directory to store the xml files
if (dir.exists("master")) {
  unlink("master", recursive = TRUE)
}
dir.create("master")

# Initialize the parameters file
param_file <- file('SIS_comp_simulations.txt', 'w')
writeLines('transmission\trecovery\tsampling\tpopulation_size', param_file)

recovery_rate <- 1

# Loop over 50 paired runs
for (i in 1:50) {
  # Randomly sample shared parameters
  # Lognormal with median ~2.5, sdlog=0.5 gives range ~1.2 to ~5.2
  transmission <- rlnorm(1, meanlog = 0.9162907318741551, sdlog = 0.5)
  while (transmission < 1.2) {
    transmission <- rlnorm(1, meanlog = 0.9162907318741551, sdlog = 0.5)
  }
  sampling <- rlnorm(1, meanlog = -4.605170185988091, sdlog = 0.25)
  population_size <- round(rlnorm(1, meanlog = 9.5, sdlog = 0.1))

  # Write parameters
  cat(sprintf('%f\t%f\t%f\t%f\n', transmission, recovery_rate, sampling, population_size), file = param_file)

  # --- Generate and run SIR simulation ---
  f <- file(sprintf('master/SIR_simulations_%d.xml', i), 'w')
  template <- file('simulation_template.xml', 'r')
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_transmission', line)) {
      writeLines(gsub('insert_transmission', as.character(transmission), line), f)
    } else if (grepl('insert_recovery', line)) {
      writeLines(gsub('insert_recovery', as.character(recovery_rate), line), f)
    } else if (grepl('insert_sampling', line)) {
      writeLines(gsub('insert_sampling', as.character(sampling), line), f)
    } else if (grepl('insert_population_size', line)) {
      writeLines(gsub('insert_population_size', as.character(population_size), line), f)
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)

  system(sprintf('/Applications/BEAST\\ 2.7.7/bin/beast -seed %d -overwrite master/SIR_simulations_%d.xml', i, i))

  # --- Generate and run SIS simulation ---
  f <- file(sprintf('master/SIS_simulations_%d.xml', i), 'w')
  template <- file('SIS_simulation_template.xml', 'r')
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_transmission', line)) {
      writeLines(gsub('insert_transmission', as.character(transmission), line), f)
    } else if (grepl('insert_recovery', line)) {
      writeLines(gsub('insert_recovery', as.character(recovery_rate), line), f)
    } else if (grepl('insert_sampling', line)) {
      writeLines(gsub('insert_sampling', as.character(sampling), line), f)
    } else if (grepl('insert_population_size', line)) {
      writeLines(gsub('insert_population_size', as.character(population_size), line), f)
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)

  system(sprintf('/Applications/BEAST\\ 2.7.7/bin/beast -seed %d -overwrite master/SIS_simulations_%d.xml', i, i))
}

close(param_file)
