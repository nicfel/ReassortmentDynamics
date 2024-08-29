library(stringr)
library(seqinr)
library(ggplot2)
library(gridExtra)

# Clear workspace
rm(list=ls())

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

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

h1n1color="#d95f02"
h3n2color="#1b9e77"


# build a wgs xml file for each virus and year between 2015 and 2018
virus = c('H1N1', 'H3N2')
year = seq(2015,2023,1)

files <- list.files("./cdcData/", pattern="*.csv", full.names=TRUE) # Get all the CSV files in cdcData
data_list <- list() # Initialize an empty list to store data frames
for (file in files) {# Read each CSV file and store it in the list
  data <- read.csv(file, header=TRUE)
  data_list[[length(data_list) + 1]] <- data
}
all_columns <- unique(unlist(lapply(data_list, colnames)))# Find all unique column names across all files
data_list <- lapply(data_list, function(df) {# Ensure each data frame has the same columns, filling in missing columns with NA
  df[setdiff(all_columns, colnames(df))] <- NA
  return(df[all_columns])
})
cases <- do.call(rbind, data_list) # Combine all data frames into one

# do the same for all csv files in CDC data/ILI
files <- list.files("./cdcData/ILI/", pattern="*.csv", full.names=TRUE)
data_list <- list()
for (file in files) {
  data <- read.csv(file, header=TRUE, skip=1)
  data_list[[length(data_list) + 1]] <- data
}
all_columns <- unique(unlist(lapply(data_list, colnames)))
data_list <- lapply(data_list, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA
  return(df[all_columns])
})
ili <- do.call(rbind, data_list)

# the cases are in the format of year, week, cases we want to convert this to a date
cases$date = as.Date(paste(cases$YEAR, cases$WEEK, 1, sep = "-"), format = "%Y-%W-%u")
cases[is.na(cases$date),]$date = as.Date("2021-01-04", format = "%Y-%m-%d")
cases = cases[order(cases$date),]

# write to file
write.csv(cases, file = "./InfA_Cases.csv", row.names = FALSE)

# the ili data is in the format of year, week, ili we want to convert this to a date
ili$date = as.Date(paste(ili$YEAR, ili$WEEK, 1, sep = "-"), format = "%Y-%W-%u")
ili[is.na(ili$date),]$date = as.Date("2021-01-04", format = "%Y-%m-%d")
ili = ili[order(ili$date),]

# go over every row and assign the cases in A..Subtyping.not.Performed., randomly
# to H3N2 or H1N1, based on the frequency of H3N2 and H1N1 in that row
for (i in seq(1, length(cases$YEAR))){
  # get the frequncy of H3N2 and H1N1 in that row
  freqH3N2 = cases$A..H3.[i]/(cases$A..H3.[i] + cases$A..2009.H1N1.[i])
  freqH1N1 = cases$A..2009.H1N1.[i]/(cases$A..H3.[i] + cases$A..2009.H1N1.[i])
  # get the number of cases in A..Subtyping.not.Performed.
  notassigned = cases$A..Subtyping.not.Performed.[i]
  # assign the cases randomly to H3N2 or H1N1
  if (notassigned == 0){
    next
  }
  if (cases$A..H3.[i]==0 & cases$A..2009.H1N1.[i]==0){
    next
  }
  addedToH3 = rbinom(1, notassigned, freqH3N2)
  cases$A..H3.[i] = cases$A..H3.[i]+ addedToH3
  cases$A..2009.H1N1.[i] = cases$A..2009.H1N1.[i]+ notassigned - addedToH3
}
# add the ili data to the cases data
cases$ili = NA
for (i in seq(1, length(cases$YEAR))){
  ilirow = ili[ili$date == cases$date[i],]
  if (nrow(ilirow) == 0){
    next
  }
  cases$ili[i] = ilirow$X..WEIGHTED.ILI
}

ili_scaler=700
p = ggplot(cases, aes(x = date)) + 
  # make a dashhe line in black with total inf A cases
  geom_line(aes(y=A..H3.+A..2009.H1N1. +B, color="Total"), linetype="dashed", color="black") +
  # geom_line(aes(y=A..H3., color="H3N2")) +
  # geom_line(aes(y=A..2009.H1N1., color="H1N1")) + 
  # add the ili data using a second axis
  geom_line(aes(x=date, y=ili*ili_scaler-(ili_scaler*1), color="ILI")) +
  theme_minimal() + 
  labs(x = "Date", y = "CDC influenza cases", title = "")
plot(p)

# treat the ili data as a function of x*totcases + y(t), estimate x using a glm
# with a non parametric y, plot x and y
# Required Libraries
library(mgcv)  # For non-parametric regression with splines
# Combine the total influenza A cases
cases$total_cases <- cases$A..H3. + cases$A..2009.H1N1. + cases$B

# Fit a Generalized Linear Model (GLM)
glm_model <- glm(ili ~ total_cases, data = cases, family = poisson())

# Extract the coefficients from the GLM model
coefficients <- coef(glm_model)

# Create a data frame to store the coefficients
coef_df <- data.frame(
  Term = names(coefficients),
  Coefficient = coefficients,
  Index = seq_along(coefficients)
)

# Plot the coefficients
library(ggplot2)

p_coefficients <- ggplot(coef_df, aes(x = Index, y = Coefficient, color = Term)) +
  geom_point(size = 3) +
  geom_line(aes(group = Term)) +
  theme_minimal() +
  labs(x = "Index", y = "Coefficient Value", title = "GLM Model Coefficients")

# Display the plot
plot(p_coefficients)

cases$predicted_ili <- predict(glm_model, type = "response")
p_comparison <- ggplot(cases, aes(x = date, y = predicted_ili)) +
  geom_line(color = "red") +
  geom_line(aes(y = ili), color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "True ILI", y = "Predicted ILI", title = "Predicted vs. True ILI")
# Display the plot
plot(p_comparison)

# compute the ratio between predicted and true ILI
cases$predicted_ili_ratio = cases$predicted_ili/cases$ili



# compute the ratio of peak total cases and peak ILI for each season
peakILIratio = data.frame()
for (y in year){
  from = as.Date(paste(y, 6, 1, sep = "-"), format = "%Y-%m-%d")
  to = as.Date(paste(y+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  peakILI = max(ili$X..WEIGHTED.ILI[ili$date>from & ili$date<to])
  peakTotal = max(cases$A..H3.[cases$date>from & cases$date<to] + cases$A..2009.H1N1.[cases$date>from & cases$date<to])
  peakILIratio = rbind(peakILIratio, data.frame(year = y, peakILI = peakILI, peakTotal = peakTotal, ratio = peakILI/peakTotal))
}
# remove 2020
peakILIratio = peakILIratio[peakILIratio$year != 2020,]

# normalize the peakILIratio to the mean value, ignoring 2020
peakILIratio$normalized = peakILIratio$ratio/mean(peakILIratio$ratio)

for (i in seq(1,length(cases$YEAR))){
  cases$cumulativeH3N2[i] = sum(cases$A..H3.[1:i])
  cases$cumulativeH1N1[i] = sum(cases$A..2009.H1N1.[1:i])
  
  # get the date of cased$date[i] and round down to the next year-05-31
  currdate = cases$date[i]
  
  #round currdate down to (year-1)-05-31 if currdate is after year-05-31
  if (currdate > as.Date(paste(cases$YEAR[i], 5, 31, sep = "-"), format = "%Y-%m-%d")){
    rounddate = as.Date(paste(cases$YEAR[i], 5, 31, sep = "-"), format = "%Y-%m-%d")
  }else{
    rounddate = as.Date(paste(cases$YEAR[i]-1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  }
  cases$cumulativeH3N2[i] = sum(cases$A..H3.[cases$date < currdate]) - sum(cases$A..H3.[cases$date < rounddate])
  cases$cumulativeH1N1[i] = sum(cases$A..2009.H1N1.[cases$date < currdate]) - sum(cases$A..2009.H1N1.[cases$date < rounddate])
  if (cases$cumulativeH3N2[i]<0){
    das
  }
}

# do the same, but read in the reassortment rates from the log files of all runs
# with constant in the name, use the column reassortmentRate.1
files <- list.files("./out2/", pattern="*.log", full.names=TRUE)
# only keep the ones that say constant
files = files[grep("constant", files)]
dynamics.constant = data.frame()
for (i in seq(1, length(files))) {
  # read in the corresponding xml file in the xmls folder that ends in xml instead of log as an xml file line by line
  filenmae = gsub("out2", "xmls", files[i])
  xmlFile = readLines(gsub("log", "xml", filenmae))
  # look for the line that <stateNode id="rateShifts" spec="RealParameter" and get the values
  rateShifts = strsplit(xmlFile[grep("<stateNode id=\"rateShifts\" spec=\"RealParameter\"", xmlFile)], 
                        split="\"")[[1]][6]
  # get the index of the line with <traitSet spec="TraitSet" traitname="date-forward" id="traitSet" dateFormat="yyyy-M-dd">
  dateline = grep("<traitSet spec=\"TraitSet\" traitname=\"date-forward\" id=\"traitSet\" dateFormat=\"yyyy-M-dd\">", xmlFile)
  # get all dates in xmlFile[dateline+1] of form yyyy-mm-dd
  datechar = xmlFile[dateline+1]
  # regexp for all yyyy-mm-dd
  dates = str_extract_all(datechar, "\\d{4}-\\d{2}-\\d{2}")[[1]]
  mrsi = max(as.Date(dates))
  # get the line with <rateShifts id="rateShifts2" spec="RealParameter" value="
  rateShifts2 = strsplit(xmlFile[grep("<rateShifts id=\"rateShifts2\" spec=\"RealParameter\" value=\"", xmlFile)], 
                         split="\"")[[1]][6]
  
  rateShiftsNe = as.numeric(strsplit(rateShifts, split=" ")[[1]])
  rateShiftsRea = as.numeric(strsplit(rateShifts2, split=" ")[[1]])
  
  
  
  t <- read.table(files[i], header = TRUE, sep = "\t")
  # take a 10% burnin
  t = t[round(nrow(t)*0.5):nrow(t),]
  # get the virus and the year from the tree file
  virus = str_extract(basename(files[i]), "H1N1|H3N2")
  year = as.numeric(str_extract(basename(files[i]), "2015|2016|2017|2018|2019|2020|2021|2022|2023"))
  
  # convert the year to a date
  date = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  rate = t$reassortmentRate.1
  HpdRate = quantile(rate, c(0.025, 0.975))
  
  # get the corresponding cummulative incidence for the same virus and year
  from = as.Date(paste(year, 6, 1, sep = "-"), format = "%Y-%m-%d")
  to = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")

  for (j in seq(1, length(rateShiftsRea)-5)){
    vals = exp(t[,paste("reassortmentRate.", j, sep="")])

    dynamics.constant = rbind(dynamics.constant, data.frame(virus = virus,
                                          time=mrsi-rateShiftsRea[j]*365-0.001,
                                          lower.5 = quantile(vals, 0.25),
                                          upper.5 = quantile(vals, 0.75),
                                          lower = quantile(vals, 0.025),
                                          upper = quantile(vals, 0.975),
                                          season=year
    ))
    dynamics.constant = rbind(dynamics.constant, data.frame(virus = virus,
                                                            time=mrsi-rateShiftsRea[j+1]*365,
                                                            lower.5 = quantile(vals, 0.25),
                                                            upper.5 = quantile(vals, 0.75),
                                                            lower = quantile(vals, 0.025),
                                                            upper = quantile(vals, 0.975),
                                                            season=year
    ))
  }
}



casesframe = data.frame(x=cases$date, y=cases$A..H3.*cases$predicted_ili_ratio, virus="H3N2")
casesframe = rbind(casesframe, data.frame(x=cases$date, y=cases$A..2009.H1N1.*cases$predicted_ili_ratio, virus="H1N1"))
# add a season to the casesframe based on if sequences where sampled between September st of year until May 31st of year+1
casesframe$season = as.numeric(format(casesframe$x, "%Y"))
casesframe$season[format(casesframe$x, "%m") < format(as.Date("2020-06-01"), "%m")] = casesframe$season[format(casesframe$x, "%m") <format(as.Date("2020-06-01"), "%m")] - 1
casesframe$season = as.factor(casesframe$season)
dynamics.constant$season = as.factor(dynamics.constant$season)
correction = 150/14
# plot the dynamics over time for both viruses
p <- ggplot(casesframe) + 
  geom_ribbon(data=dynamics.constant, aes(x = time, ymin = lower.5, ymax = upper.5, group=season, fill=virus)) +
  geom_ribbon(data=dynamics.constant, aes(x = time, ymin = lower, ymax = upper, group=season, fill=virus), alpha = 0.5, color = NA) +
  # Add the cases over time
  geom_line(aes(x = x, y = y/8000), color="black", linetype="dashed") +  # Fixed column reference with backticks
  facet_grid(virus~season, scales = "free_x") +  # Fixed the facet_wrap formula
  coord_cartesian(ylim = c(0, 0.5)) +  # Fixed the coord_cartesian formula
  theme_minimal() + 
  scale_fill_manual(values=c("H1N1"=h1n1color, "H3N2"=h3n2color))+
  labs(x = "Date", y = "Reassortment Rate", title = "")
# Plot the graph
print(p)
ggsave("./../../Figures/InfA_Prevalence_piecewise.pdf", p, width = 9, height = 5)

# do the same, but read in the reassortment rates from the log files of all runs
# with constant in the name, use the column reassortmentRate.1
files <- list.files("./out2/", pattern="*.log", full.names=TRUE)
# only keep the ones that say constant
files = files[grep("variable", files)]
dynamics = data.frame()
ne = data.frame()
for (i in seq(1, length(files))) {
  # read in the corresponding xml file in the xmls folder that ends in xml instead of log as an xml file line by line
  filenmae = gsub("out2", "xmls", files[i])
  xmlFile = readLines(gsub("log", "xml", filenmae))
  # look for the line that <stateNode id="rateShifts" spec="RealParameter" and get the values
  rateShifts = strsplit(xmlFile[grep("<stateNode id=\"rateShifts\" spec=\"RealParameter\"", xmlFile)], 
                                   split="\"")[[1]][6]
  # get the index of the line with <traitSet spec="TraitSet" traitname="date-forward" id="traitSet" dateFormat="yyyy-M-dd">
  dateline = grep("<traitSet spec=\"TraitSet\" traitname=\"date-forward\" id=\"traitSet\" dateFormat=\"yyyy-M-dd\">", xmlFile)
  # get all dates in xmlFile[dateline+1] of form yyyy-mm-dd
  datechar = xmlFile[dateline+1]
  # regexp for all yyyy-mm-dd
  dates = str_extract_all(datechar, "\\d{4}-\\d{2}-\\d{2}")[[1]]
  mrsi = max(as.Date(dates))
  # get the line with <rateShifts id="rateShifts2" spec="RealParameter" value="
  rateShifts2 = strsplit(xmlFile[grep("<rateShifts id=\"rateShifts2\" spec=\"RealParameter\" value=\"", xmlFile)], 
                        split="\"")[[1]][6]
  
  rateShiftsNe = as.numeric(strsplit(rateShifts, split=" ")[[1]])
  rateShiftsRea = as.numeric(strsplit(rateShifts2, split=" ")[[1]])
  
  
  
  t <- read.table(files[i], header = TRUE, sep = "\t")
  # take a 10% burnin
  t = t[round(nrow(t)*0.7):nrow(t),]
  # get the virus and the year from the tree file
  virus = str_extract(basename(files[i]), "H1N1|H3N2")
  year = as.numeric(str_extract(basename(files[i]), "2015|2016|2017|2018|2019|2020|2021|2022|2023"))

  # convert the year to a date
  date = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  rate = t$reassortmentRate.1
  HpdRate = quantile(rate, c(0.025, 0.975))
  
  # get the corresponding cummulative incidence for the same virus and year
  from = as.Date(paste(year, 6, 1, sep = "-"), format = "%Y-%m-%d")
  to = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")

  timestep = max(rateShiftsNe)/1000
  no_points = floor(0.7/timestep)
  
  # dynamics = rbind(dynamics, data.frame(virus = virus,
  #                                       time=mrsi-(timestep)/2*365,
  #                                       median = 0, lower = 0, upper = 0
  # ))
  
  
  for (j in seq(0, no_points)){
    vals = t[,paste("reassortment", j, sep="")]
    dynamics = rbind(dynamics, data.frame(virus = virus,
                                          time=mrsi-timestep*j*365,
                                          lower.5 = quantile(vals, 0.25),
                                          upper.5 = quantile(vals, 0.75),
                                          lower = quantile(vals, 0.025),
                                          upper = quantile(vals, 0.975),
                                          season=year
                                          ))
  }
  
  
  for (i in 1:length(rateShiftsNe)) {
    # start a new matrix of size length(unique(t$Sample)) and 1000
    I = matrix(0, nrow=length(t$Sample), ncol=1000)
    # init a length(timepoints)x4 matrix for the splineCoefficents
    splineCoeffs = matrix(0, nrow=length(rateShiftsNe)-1, ncol=4)
    # for each Sample in t, compute the Ne trajectory
    for (s in 1:length(t$Sample)) {
      # populate the spline coeffients for this iteration of t of the names splineCoeffs_0_0....
      for (a in 1:length(rateShiftsNe)-1){
        for (b in seq(1, 4)){
          splineCoeffs[a,b] = t[s, paste("splineCoeffs", a-1, b-1, sep="_")]
        }
      }
      I[s,] = interpolate_I_over_grid(rateShiftsNe, splineCoeffs, 1000, max(rateShiftsNe)/1000)$I
    }
    # loop over all colums in I
    for (j in seq(1,no_points,5)) {
      # add the time and rate to the ne dataframe
      ne = rbind(ne, data.frame(virus = virus,
                                time=mrsi-j*max(rateShiftsNe)/1000*365,
                                upper.5 = quantile(I[,j], 0.75),
                                lower.5 = quantile(I[,j], 0.25),
                                upper = quantile(I[,j], 0.975),
                                lower = quantile(I[,j], 0.025),
                                season=year))
    }
  }
  
  
  
  
  
  # add 0 values after
  # dynamics = rbind(dynamics, data.frame(virus = virus,
  #                                       time=mrsi-timestep*(no_points+1)*365,
  #                                       median = 0, lower = 0, upper = 0
  # ))
}


correction = 150/14
# plot the dynamics over time for both viruses
p <- ggplot(casesframe) + 
  geom_ribbon(data=dynamics, aes(x = time, ymin = lower.5, ymax = upper.5, group=season, fill=virus)) +
  geom_ribbon(data=dynamics, aes(x = time, ymin = lower, ymax = upper, group=season, fill=virus), alpha = 0.5, color = NA) +
  # Add the cases over time
  geom_line(aes(x = x, y = y/4000), color="black", linetype="dashed") +  # Fixed column reference with backticks
  facet_grid(virus~.) +  # Fixed the facet_wrap formula
  coord_cartesian(ylim = c(0, 1)) +  # Fixed the coord_cartesian formula
  theme_minimal() + 
  labs(x = "Date", y = "Reassortment Rate", title = "")+
  scale_fill_manual(values=c("H1N1"=h1n1color, "H3N2"=h3n2color))
# Plot the graph
print(p)
ggsave("./../../Figures/InfA_Prevalence.pdf", p, width = 9, height = 5)


# plot the dynamics over time for both viruses
p <- ggplot(casesframe) + 
  geom_ribbon(data=ne, aes(x = time, ymin = lower.5, ymax = upper.5, group=season, fill=virus)) +
  geom_ribbon(data=ne, aes(x = time, ymin = lower, ymax = upper, group=season, fill=virus), alpha = 0.5, color = NA) +
  # Add the cases over time
  geom_line(aes(x = x, y = y/5), color="black", linetype="dashed") +  # Fixed column reference with backticks
  facet_grid(virus~.) +  # Fixed the facet_wrap formula
  coord_cartesian(ylim = c(0, 1000)) +  # Fixed the coord_cartesian formula
  theme_minimal() + 
  labs(x = "Date", y = "Ne", title = "")+
  scale_fill_manual(values=c("H1N1"=h1n1color, "H3N2"=h3n2color))
# Plot the graph
print(p)
ggsave("./../../Figures/InfA_Ne.pdf", p, width = 9, height = 5)



