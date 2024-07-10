# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# Define superinfection exclusion dynamics
percentage_success <- c(0.5, 0.5, 0.5, 0.3, 0.25, 0.25, 0.1, 0.025, 0, 0)
time <- c(0:(length(percentage_success) - 2), 10000)

# Define the recovery rate as an exponential decay with mean 4 days.
recovery_rate <- 1 / 5

# Assume a constant rate of infection
infection_rate <- 1

# Define the number of simulations
n_simulations <- 1000000
totProb <- 0

# Run simulations
for (r in 1:n_simulations) {
  # Define the simulation time as an exponential random variable with rate recovery_rate 
  simulationTime <- rexp(1, rate = recovery_rate) * 24
  
  # Define what percentage of infections would have led to a superinfection before the first infection is cleared
  totPercentage <- 0
  i <- 1
  while (time[i + 1] < simulationTime) {
    totPercentage <- totPercentage + percentage_success[i] * (time[i + 1] - time[i])
    i <- i + 1
  }
  totPercentage <- totPercentage + percentage_success[i] * (simulationTime - time[i])
  totProb <- totProb + totPercentage / simulationTime
}

SIE_prob <- 1 / (totProb / n_simulations)

# Generate samples and prevalence values in log space 
nsamples <- 10^(seq(1, 5, length.out = 100))
prevalence <- 10^(seq(-5, -1, length.out = 100))
plot_data = data.frame()
# compute meanEvents <- nsamples * prevalence * 150 / 14 * 7 / 365)
# fore every combination of nsamples and prevalence
for (i in 1:length(nsamples)) {
  for (j in 1:length(prevalence)) {
    meanEvents <- nsamples[i] * prevalence[j] * 150 / 14 * 7 / 365
    plot_data <- rbind(plot_data, data.frame(nsamples = nsamples[i], prevalence = prevalence[j], meanEvents = meanEvents))
  }
}

# compute lines for meanEvents = 1, 10, 100, 1000
lines_data = data.frame()
targets = c(1, 10, 100, 1000)
for (i in 1:length(targets)) {
  for (j in 1:length(nsamples)) {
    prev <- targets[i] / (nsamples[j] * 150 / 14 * 7 / 365)
    lines_data <- rbind(lines_data, data.frame(nsamples = nsamples[j], prevalence = prev, target=targets[i]))
  }
}

# remove any values outside max(prevalence) or max(nsamples)
lines_data <- lines_data[lines_data$prevalence <= max(prevalence),]

# make a dataframe to label the lines using the mid point for nsample and prevalence for each target
text_dat = data.frame()
for (t in targets){
  tmp = lines_data[lines_data$target == t,]
  text_dat = rbind(text_dat, data.frame(nsamples = median(tmp$nsamples), prevalence = median(tmp$prevalence), target = t))
}
                        
# plot the nsamples on the x, the prevalence on the y and color each tile by meanEvents
p = ggplot(plot_data, aes(x=nsamples, y=prevalence, fill=log(meanEvents)))+
  geom_tile()+
  scale_fill_viridis_c()+
  scale_x_log10()+
  scale_y_log10()+
  labs(x = "Number of samples", y = "Prevalence", fill = "log number of expected\nreassortment events)")+
  # add lines for meanEvents = 1, 10, 100, 1000
  geom_line(data = lines_data, aes(x = nsamples, y = prevalence, group=target), color = "grey")+
  geom_text(data = text_dat, aes(x = nsamples, y = prevalence, label = target), color = "black")+
  theme_minimal()
plot(p)
