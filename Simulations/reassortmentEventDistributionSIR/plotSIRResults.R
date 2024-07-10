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

# set the time points for the parameters
time_points = seq(0,12,1.5)
# get 500 points between 0 and 16
time_points2 = seq(0,12, length.out=500)

# read in the file SIR_simulations.txt that contains the simulated values
simulated = read.table("SIR_simulations.txt", header=TRUE, sep="\t")


# get the true values for the prevalence by readin in the files in /master
true.data = data.frame()
# get all .log files in .master
logfiles <- list.files("./master", pattern="*.log", full.names=TRUE)
# loop over all files
for (i in seq(1, length(logfiles))){
# for (i in seq(1,10)){
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][4])
  # read in the file line by line
  lines = readLines(logfiles[i])
  data = strsplit(lines[[2]], split="\t")[[1]]

  print(logfiles[i])
  
  # Initialize lists to hold the results
  I_list <- list()
  time_I_list <- list()
  co_inf_times_list <- list()

  # Replace [ and ] then split the SIR events into a vector of strings
  sir <- str_split(str_replace_all(data[[2]], "\\[|\\]", ""), ', ')[[1]]

  # Compute the number of infected individuals over time
  I <- 1
  time_I <- 0
  # Keep also track of the co-infection event times
  total_infected = 0
  total_coinfection = 0

  for (j in seq_along(sir)) {
      tmp <- str_split(sir[j], ':')[[1]]
      if (tmp[1] == '0') {
          I <- c(I, I[length(I)])
          # time_I <- c(time_I, as.numeric(tmp[2]))
          total_coinfection = total_coinfection+1
      } else if (tmp[1] == '1') {
          total_infected = total_infected+1
          I <- c(I, I[length(I)] + 1)
          # time_I <- c(time_I, as.numeric(tmp[2]))
      } else if (tmp[1] == '2') {
          I <- c(I, I[length(I)] - 1)
          # time_I <- c(time_I, as.numeric(tmp[2]))
      } else if (tmp[1] == '3') {
          I <- c(I, I[length(I)])
          # time_I <- c(time_I, as.numeric(tmp[2]))
      }
  }
  
  # get the population size
  popSize = simulated[run, "population_size"]
  recovered = 
  # add the data to the data frame
  true.data = rbind(true.data, data.frame(run=run, avg_prev = mean(I),
                                          totalInfections=total_infected, 
                                          transmissionRate = simulated[run, "transmission"],
                                          popSize=popSize, 
                                          recovered=simulated[run, "recovery"],
                                          total_coinfection=total_coinfection))
}


p1=ggplot(true.data, aes(x=recovered, y=total_coinfection/popSize))+
  geom_point()+
  geom_smooth()
p2=ggplot(true.data, aes(x=totalInfections/popSize, y=total_coinfection/popSize))+
  geom_point()+
  geom_smooth()
p3=ggplot(true.data, aes(x=avg_prev/popSize, y=total_coinfection/popSize))+
  geom_point()

ggplot(true.data, aes(x=recovered, y=transmissionRate,z=total_coinfection/totalInfections))+
  geom_point(aes(color=log(total_coinfection/avg_prev)))
  # geom_contour_filled(aes(color=..level..))

p = gridExtra::grid.arrange(p1, p2, p3, ncol=3)
plot(p)



