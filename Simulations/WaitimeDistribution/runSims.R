# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# to keep me from accidentally running these
# dasda
# dasdas
# das
# das

# clear workspace
rm(list = ls())

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
writeLines('transmission\trecovery\tsampling\tpopulation_size\tk', param_file)

recovery_rate <- 1

# Loop over 100 runs
for (i in 1:1000) {
  # Make a new xml file for each run
  f <- file(sprintf('master/SIR_simulations_%d.xml', i), 'w')
  # Open the template file
  template <- file('../SIR/simulation_template.xml', 'r')
  
  # Randomly sample the transmission rate from a lognormal distribution with mean 3 and S 0.25
  transmission <- rlnorm(1, meanlog = 1.0986122886681098, sdlog = 0.25)
  while (transmission < 1.2) {
    transmission <- rlnorm(1, meanlog = 1.0986122886681098, sdlog = 0.25)
  }
  
  # Randomly sample the sampling rate from a lognormal distribution with mean 0.01 and S 0.5
  sampling <- rlnorm(1, meanlog = -4.605170185988091, sdlog = 0.25)
  # Randomly sample the population size from a lognormal distribution with mean 10000 and S 0.5
  population_size <- round(rlnorm(1, meanlog = 9.5, sdlog = 0.1))
  # Randomly sample the k of the negative binomial distribution from a lognormal distribution with mean 1 and S 0.5
  k <- rlnorm(1, meanlog = 0, sdlog = 1)

  
  # Write the parameters to the file
  cat(sprintf('%f\t%f\t%f\t%f\t%f\n', transmission, recovery_rate, sampling, population_size, k), file=param_file)
  
  # Write the parameters to the xml file
  if (i>500){
    while (length(line <- readLines(template, n = 1)) > 0) {
      if (grepl('spec="SIRwithReassortment"', line)) {
        writeLines(gsub('spec="SIRwithReassortment"', 'spec="SuperspreadingSIRwithReassortment"', line), f)
      } else if(grepl('insert_transmission', line)) {
        writeLines(gsub('insert_transmission', as.character(transmission), line), f)
      } else if (grepl('insert_recovery', line)) {
        writeLines(gsub('insert_recovery', as.character(recovery_rate), line), f)
      } else if (grepl('insert_sampling', line)) {
        writeLines(gsub('insert_sampling', as.character(sampling), line), f)
      } else if (grepl('insert_population_size', line)) {
        writeLines(gsub('insert_population_size', as.character(population_size), line), f)
        writeLines(gsub('populationSize="insert_population_size"', paste('k="', as.character(k), '"', sep=""), line), f)
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)
  }else{
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
  }
  
  # Run the xml using BEAST and the system command
  system(sprintf('/Applications/BEAST\\ 2.7.6/bin/beast -seed %d -overwrite master/SIR_simulations_%d.xml', i, i))
}

close(param_file)

