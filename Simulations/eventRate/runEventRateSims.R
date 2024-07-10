library(stringr)
# Clear workspace
rm(list=ls())

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

# Initialize the parameters for the SIR simulations
param_file <- file('SIR_simulations.txt', 'w')

# define the number of states
nstates = 20

# Initialize the parameters for the SIR simulations
param_file <- file('eventRateSIR_simulations.txt', 'w')

# write one header with a different rate for each state from 1...nstates

writeLines("transmission\trecovery\tpopulation_size\td", param_file)
recovery_rate <- 1

# Loop over 100 runs
for (i in 1:100) {
  print(i)
  # Make a new xml file for each run
  f <- file(sprintf('master/eventRateSIR_simulations_%d.xml', i), 'w')
  # Open the template file
  template <- file('../structuredSIR/structuredSimulation_template.xml', 'r')

  # Randomly sample one transmission rate from a lognormal distribution with mean in real space of 1.5 and S 0.25
  transmission <- rlnorm(1, meanlog = 0.4, sdlog = 0.25)
  while (transmission < 1.01) {
    transmission <- rlnorm(1, meanlog = 0.4, sdlog = 0.25)
  }

  # Randomly sample 10 sampling rates from a lognormal distribution with mean in real space of 0.01 and S 0.5
  sampling <- rlnorm(nstates, meanlog = -3.605170185988091, sdlog = 0.5)
  # Randomly sample 10 population sizes from a lognormal distribution with mean 10000 and S 0.5
  population_size <- round(rlnorm(nstates, meanlog = 6.510340371976182, sdlog = 0.5))
  # Randomly sample the k of the negative binomial distribution from a lognormal distribution with mean 1 and S 0.5
  k <- rlnorm(1, meanlog = 0, sdlog = 1)
  # Randomly sample nstates*(nstates-1) migration rates from a log normal with mean in real space of 0.1 and S of 0.25
  if (i>50){
    migration <- rlnorm(nstates*(nstates-1), meanlog = -2.5, sdlog = 1)
  }else{
    migration <- rlnorm(nstates*(nstates-1), meanlog = -4.5, sdlog = 1)
  }

  # replace transmission, recover and k and waning (just 0's) nstates times
  transmission = rep(transmission, nstates)
  recovery = rep(recovery_rate, nstates)
  waning = rep(0,nstates)

  # Write the parameters to the file
  cat(sprintf('%f\t%f\t%f\t%f\n', transmission[[1]], recovery_rate[[1]], sum(population_size), k), file=param_file)

  # Write the parameters to the xml file
  while (length(line <- readLines(template, n = 1)) > 0) {
    if(grepl('insert_transmission', line)) {
      writeLines(gsub('insert_transmission', paste(transmission, collapse=" "), line), f)
    } else if (grepl('insert_recovery', line)) {
      writeLines(gsub('insert_recovery', paste(recovery, collapse=" "), line), f)
    } else if (grepl('insert_sampling', line)) {
      writeLines(gsub('insert_sampling', paste(sampling, collapse=" "), line), f)
    } else if (grepl('insert_population_size', line)) {
      writeLines(gsub('insert_population_size', paste(population_size, collapse=" "), line), f)
    } else if (grepl('insert_migration', line)) {
      writeLines(gsub('insert_migration', paste(migration, collapse=" "), line), f)
    } else if (grepl('insert_waning', line)) {
      writeLines(gsub('insert_waning', paste(waning, collapse=" "), line), f)
    }else if (grepl('insert_k', line)){
      writeLines(gsub('insert_k', as.character(k), line), f)
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)
  # Run the xml using BEAST and the system command, while preventing any logging to screen
  system(sprintf('/Applications/BEAST\\ 2.7.6/bin/beast -seed %d -overwrite master/eventRateSIR_simulations_%d.xml', i, i))
}
close(param_file)

