library(stringr)
library(seqinr)
library(ggplot2)
library(cowplot)


# 2nd part of plotDynamics.R
segments = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

# Set the directory to the directory of the file
setwd("~/wildbird_HPAI_LPAI/")
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# function that redoes the spline interpolation that is directly done in CoalRe
interpolate_I_over_grid <- function(rateShifts, splineCoeffs, gridPoints, dt) {
  # Initialize variables
  time <- numeric(gridPoints)
  I <- numeric(gridPoints)
  j <- 1
  k <- j-1
  
  # Loop over grid points
  for (i in 1:gridPoints) {
    # Update the time for this grid point
    time[i] <- (i - 1) * dt
    
    
    # Find the interval in which this grid point lies
    if (time[i] >= rateShifts[j]) {
      j <- j + 1
      k <- k + 1
      if (k == length(rateShifts)) {
        k <- k - 1
      }
    }
    
    # Get the time diff from the last point where logI was estimated
    
    timeDiff <- time[i] - rateShifts[k]
    timeDiff2 <- timeDiff^2
    timeDiff3 <- timeDiff2 * timeDiff
    if (timeDiff<0){
      shouldnothappen
    }
    
    # Compute the number of infected individuals at the grid points
    I[i] <- exp(splineCoeffs[k, 1] * timeDiff3 + splineCoeffs[k, 2] * timeDiff2 + splineCoeffs[k, 3] * timeDiff + splineCoeffs[k, 4])
  }
  return(list(time = time, I = I))
}



data = data.frame()

burnin = 0.8

# get all files TODO: the file names should be replace with the variable ones
t1 = read.table("variable/rep3/H5N1_wgs_HPAI_LPAI_North_America.variable.rep0.log", header = TRUE, sep = "\t")
t2 = read.table("variable/rep2/H5N1_wgs_HPAI_LPAI_North_America.variable.rep2.log", header = TRUE, sep = "\t")
t3 = read.table("variable/rep1/H5N1_wgs_HPAI_LPAI_North_America.variable.rep1.log", header = TRUE, sep = "\t")

# combined after 10% burnin
t = rbind(t1[round(nrow(t1)*burnin):nrow(t1),], t2[round(nrow(t2)*burnin):nrow(t2),], t3[round(nrow(t3)*burnin):nrow(t3),])
# read in the xml file xmls/H5N1_wgs.varying.rep0.xml
# TODO change to the corresponding xml file
xmlFile = readLines("variable/rep1/H5N1_wgs_HPAI_LPAI_North_America.variable.rep1.xml")
# look for the line that <stateNode id="rateShifts" spec="RealParameter" and get the values
rateShifts = strsplit(xmlFile[grep("<stateNode id=\"rateShifts\" spec=\"RealParameter\"", xmlFile)], 
                      split="\"")[[1]][6]
# get the index of the line with <traitSet spec="TraitSet" traitname="date-forward" id="traitSet" dateFormat="yyyy-M-dd">
dateline = grep("<traitSet spec=\"TraitSet\" traitname=\"date-forward\" id=\"traitSet\" dateFormat=\"yyyy-M-dd\">", xmlFile)
# get all dates in xmlFile[dateline+1] of form yyyy-mm-dd
datechar = xmlFile[dateline+1]
# regexp for all yyyy-mm-dd
dates = str_extract_all(datechar, "\\d{4}-\\d{2}-\\d{2}")[[1]]
#dates = dates[!is.na(dates)]

mrsi = max(as.Date(dates))
# get the line with <rateShifts id="rateShifts2" spec="RealParameter" value="
rateShifts2 = strsplit(xmlFile[grep("<rateShifts id=\"rateShifts2\" spec=\"RealParameter\" value=\"", xmlFile)], 
                       split="\"")[[1]][6]
# clear xmlFIle from memory
rm(xmlFile)
# make two datframes, the first reassortment, reads in the piecewise constant reassortment rates using rateShifts2
reassortment = data.frame()
timepoints = as.numeric(strsplit(rateShifts2, split=" ")[[1]])
timepoints_fine = seq(min(timepoints), max(timepoints), length.out = 999)

for (i in 1:(length(timepoints_fine))) {
  rate = t[, paste("reassortment", i, sep="")]
  reassortment = rbind(reassortment, data.frame(from=mrsi-timepoints_fine[i]*365, 
                                                rate=median(rate), 
                                                upper.5=quantile(rate, 0.75), lower.5=quantile(rate, 0.25),
                                                upper=quantile(rate, 0.975), lower=quantile(rate, 0.025)))
  reassortment = rbind(reassortment, data.frame(from=mrsi-timepoints_fine[i+1]*364, 
                                                rate=median(rate), 
                                                upper.5=quantile(rate, 0.75), lower.5=quantile(rate, 0.25),
                                                upper=quantile(rate, 0.975), lower=quantile(rate, 0.025)))
  
}

# read in the effective population sizes
ne = data.frame()
timepoints = as.numeric(strsplit(rateShifts, split=" ")[[1]])
for (i in 1:length(timepoints)) {
  # start a new matrix of size length(unique(t$Sample)) and 1000
  I = matrix(0, nrow=length(t$Sample), ncol=1000)
  # init a length(timepoints)x4 matrix for the splineCoefficents
  splineCoeffs = matrix(0, nrow=length(timepoints)-1, ncol=4)
  # for each Sample in t, compute the Ne trajectory
  for (s in 1:length(t$Sample)) {
    # populate the spline coeffients for this iteration of t of the names splineCoeffs_0_0....
    for (a in 1:length(timepoints)-1){
      for (b in seq(1, 4)){
        splineCoeffs[a,b] = t[s, paste("splineCoeffs", a-1, b-1, sep="_")]
      }
    }
    I[s,] = interpolate_I_over_grid(timepoints, splineCoeffs, 1000, max(timepoints)/1000)$I
  }
  
  # loop over all colums in I
  for (j in 1:300) {
    # add the time and rate to the ne dataframe
    ne = rbind(ne, data.frame(time=mrsi-j*max(timepoints)/1000*365,
                              upper.5 = quantile(I[,j], 0.75),
                              lower.5 = quantile(I[,j], 0.25),
                              upper = quantile(I[,j], 0.975),
                              lower = quantile(I[,j], 0.025)))
  }
}


# plot the two
p1 = ggplot(ne, aes(x=time)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) + 
  geom_ribbon(aes(ymin=lower.5, ymax=upper.5), alpha=0.75) + 
  coord_cartesian(xlim=c(mrsi-365*2.5, mrsi), ylim=c(0,500)) +
  xlab("") + ylab("Effective population\nsize")+
  theme_minimal() +
  theme(text = element_text(size = 18)) 

plot(p1)


p2 = ggplot(reassortment, aes(x=from)) +
  geom_ribbon(aes(x=from, ymin=lower.5, ymax=upper.5), alpha=0.5, fill="#ff7f00") + 
  geom_ribbon(aes(x=from, ymin=lower, ymax=upper), alpha=0.75, fill="#ff7f00") + 
  ylim(0,1) +
  coord_cartesian(xlim=c(mrsi-365*2.5, mrsi)) +
  xlab("") + ylab("Reassortment rate")+
  theme_minimal() +
  theme(text = element_text(size = 18)) 


plot(p2)

# read in HPAI detection data from USDA
HPAI_pos_day_date <- read.csv("~/Documents/HPAI/HPAI_epi_casevis/all_hpai_full.csv")
HPAI_pos_day_date$Date <- as.Date(HPAI_pos_day_date$Date)

p3 = ggplot(HPAI_pos_day_date, aes(x = Date, y = count)) +
  geom_bar(stat = "identity") +
  xlab("Date") +
  ylab("Detections\nof HPAI") +
  coord_cartesian(xlim=c(mrsi-365*2.5, mrsi)) +
  theme_minimal() +
  theme(text = element_text(size = 18)) 

p3


pc <- plot_grid(p1,p2,p3, nrow = 3, align = "v")
pc


ggsave("Figures/H5N1_HPAI_LPAI_variable_dynamics.pdf", pc, width = 11, height = 8.5, units = "in")



#####
# export the dataframes

write.csv(ne,"HPAI_LPAI_variable_Ne.csv")
write.csv(reassortment,"HPAI_LPAI_variable_reassortmentrate.csv")
