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
simulated = read.table("structuredSIR_simulations.txt", header=TRUE, sep="\t")

# get all .log files in /out
logfiles <- list.files("./out", pattern="*.log", full.names=TRUE)
# loop over all files
for (i in seq(1, length(logfiles))){
  
  # split on . and get the second group as the method used
  method = strsplit(logfiles[i], "\\.")[[1]][3]
  print(method)
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][4])
  # if the method is not constant
  if (grepl("constant", method)|| grepl("variable", method)){
    next
  }
  
  # try to open logfiles[[i]] as t = read.table(logfiles[i], header=TRUE, sep="\t"), otherwise skip it
  t = try(read.table(logfiles[i], header=TRUE, sep="\t"))
  # if t is not a data.frame, skip the file
  if (!is.data.frame(t)){
    next
  }

  # remove 10 % as burnin
  t = t[round(0.1*nrow(t)):nrow(t),]
  # calculate the ESS for columns 2 to 5
  ess = sapply(t[,2:5], function(x) effectiveSize(x))
  # calculate the minimum ESS
  minESS = min(ess)
  # if the minimum ESS is below 100, print the filename
  if (minESS < 50) {
    print(minESS)
    print(logfiles[i])
    next
  }
  # get the value in simulated in the run row and "transmission" column
  transmision0 = simulated[run, "transmission"]
  
  reassortment.data = rbind(reassortment.data, data.frame(method=method, run=run, 
          transmission = simulated[run, "transmission"],
          k = simulated[run, "d"],
          events=median(t[, "InfectedToRho"]), 
          eventsl=quantile(t[, "InfectedToRho"], 0.025), 
          eventsu=quantile(t[, "InfectedToRho"], 0.975)))

}


reassortment.data = reassortment.data[reassortment.data$events>0.00000001,]
p = ggplot(reassortment.data, aes(x=transmission*(1+1/k), y=events, ymin=eventsl, ymax=eventsu))+
  geom_point()+geom_errorbar()+scale_y_log10() + scale_x_log10() + 
  ylab("constant for Ne=c*I ")+
  ylab("constant for Ne=c*I ")+ 
  theme_minimal()
  # geom_smooth(formula = "y~x+1")
ggsave("./../../Figures/StructuredSir_neRalation.pdf", p, width = 6, height = 3)


