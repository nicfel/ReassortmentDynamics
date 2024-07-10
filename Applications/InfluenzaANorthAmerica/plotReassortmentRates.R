library(stringr)
library(seqinr)
library(ggplot2)
library(ggpubr)

# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# build a wgs xml file for each virus and year between 2015 and 2018
virus = c('H1N1', 'H3N2')
year = seq(2015,2023,1)
h1n1color="#d95f02"
h3n2color="#1b9e77"

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

# the ili data is in the format of year, week, ili we want to convert this to a date
ili$date = as.Date(paste(ili$YEAR, ili$WEEK, 1, sep = "-"), format = "%Y-%W-%u")
ili[is.na(ili$date),]$date = as.Date("2021-01-04", format = "%Y-%m-%d")
ili = ili[order(ili$date),]


# plot the TOTAL.SPECIMENS over time
p = ggplot(cases, aes(x = date)) + 
  # geom_line(aes(y=TOTAL.SPECIMENS, color="Total"), linetype="dashed", color="black") +
  # geom_line(aes(y=A..H3.+A..2009.H1N1.+A..Subtyping.not.Performed.+B+H3N2v+BVic+BYam))+
  geom_line(data=ili, aes(x=date, y=X..WEIGHTED.ILI*500, color="ILI")) +
  theme_minimal() + 
  labs(x = "Date", y = "CDC influenza cases", title = "")
plot(p)


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



# read in the Global WHO influenza cases downloaded from FluNet
who_all = read.csv("./whoData/North_SW_Europe.csv", header = TRUE)
who_all$date = as.Date(who_all$ISO_SDATE)
who_NA = read.csv("./whoData/NA_who_sentinel.csv", header = TRUE)
who_NA$date = as.Date(who_NA$ISO_SDATE)
# make a new data frame who, that summarize the cases over all locations, for 
# weeks with the same date
who = data.frame()
dates = unique(who_all$date)
for (i in seq(1, length(dates))){
  h3 = who_all$AH3[who_all$date == dates[i]]
  h1 = who_all$AH1N12009[who_all$date == dates[i]]
  # sum over all non NA values
  who = rbind(who, data.frame(date = dates[i], cases=sum(h3[!is.na(h3)]), virus="H3N2", location="Western Europe"))
  who = rbind(who, data.frame(date = dates[i], cases=sum(h1[!is.na(h1)]), virus="H1N1", location="Western Europe"))
  h3 = who_NA$AH3[who_NA$date == dates[i]]
  h1 = who_NA$AH1N12009[who_NA$date == dates[i]]
  # sum over all non NA values
  who = rbind(who, data.frame(date = dates[i], cases=sum(h3[!is.na(h3)]), virus="H3N2", location="North America"))
  who = rbind(who, data.frame(date = dates[i], cases=sum(h1[!is.na(h1)]), virus="H1N1", location="North America"))
}




p = ggplot(cases, aes(x = date)) + 
  # make a dashhe line in black with total inf A cases
  geom_line(aes(y=A..H3.+A..2009.H1N1. +B, color="Total"), linetype="dashed", color="black") +
  geom_line(aes(y=A..H3., color="H3N2")) +
  geom_line(aes(y=A..2009.H1N1., color="H1N1")) +
  # add the ili data using a second axis
  geom_line(data=ili, aes(x=date, y=X..WEIGHTED.ILI*500, color="ILI")) +
  theme_minimal() + 
  labs(x = "Date", y = "CDC influenza cases", title = "")
plot(p)

p = ggplot(who, aes(x = date)) + 
  # make a dashhe line in black with total inf A cases
  geom_line(aes(y=cases, color=virus, linetype = location)) +
  theme_minimal() + 
  labs(x = "Date", y = "WHO sentinel cases", title = "")
plot(p)


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

p = ggplot(cases, aes(x = date)) + 
  geom_line(aes(y=cumulativeH3N2, color="H3N2")) + 
  geom_line(aes(y=cumulativeH1N1, color="H1N1")) + 
  theme_minimal() + 
  labs(x = "Date", y = "cumulative CDC influenza incidence until May 31st", title = "")
plot(p)

# loop over all trees files that have constant in the name
tree.files <- list.files("./out/", pattern="*.trees", full.names=TRUE)
# only keep the constant files
# tree.files = tree.files[grep("constant", tree.files)]
maxTipDistance = 12/12
for (i in seq(1, length(tree.files))) {
  print(tree.files[i])
  system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
         paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ExtractReassortmentEvents",
                         "-burnin 60 ", tree.files[i],
                         gsub(".trees", ".events.reassortment", tree.files[i])))
  
  # summarize the networks
  # system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
  #        paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentNetworkSummarizer",
  #                        "-burnin 20", tree.files[i],
  #                        gsub(".trees", ".tree", tree.files[i])))
  # 
  system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
         paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ExtantReassortmentRate",
                         "-burnin 60 -maxTipDistance", maxTipDistance, tree.files[i],
                         gsub(".trees", ".tip.reassortment", tree.files[i])))
  # for (a in seq(0,7)){
  #   system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T, paste("/Applications/BEAST\\ 2.7.6/bin/applauncher TipDistanceReassortmentRate -removeSegments",a,
  #                          "-burnin 60 -maxTipDistance", maxTipDistance, tree.files[i],
  #                          gsub(".trees", paste(".",a,".tip.reassortment", sep=""), tree.files[i])))
  # }
}

# read in all the tip.reassortment files
files <- list.files("./out/", pattern="*.tip.reassortment", full.names=TRUE)
adjusted_reassortment = data.frame()
for (i in seq(1, length(files))) {
  t <- read.table(files[i], header = FALSE, sep = "\t")

  # get the virus and the year from the tree file
  virus = str_extract(basename(files[i]), "H1N1|H3N2")
  year = as.numeric(str_extract(basename(files[i]), "2015|2016|2017|2018|2019|2020|2021|2022|2023"))
  method = str_extract(basename(files[i]), "constant|ne|variable")
  # split the file names on the dot
  tmp = strsplit(basename(files[i]), "\\.")[[1]]
  if (length(tmp) == 5){
    removedSegments = tmp[4]
  }else{
    removedSegments = "none"
  }


  # convert the year to a date
  date = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  rate = t$V3
  HpdRate = quantile(rate, c(0.025, 0.975))

  # get the corresponding cummulative incidence for the same virus and year
  from = as.Date(paste(year, 6, 1, sep = "-"), format = "%Y-%m-%d")
  to = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  if (virus == "H1N1"){
    cumulative = max(cases$cumulativeH1N1[cases$date>from & cases$date<to])
    peakCases = max(cases$A..2009.H1N1.[cases$date>from & cases$date<to])
  }else{
    cumulative = max(cases$cumulativeH3N2[cases$date>from & cases$date<to])
    peakCases = max(cases$A..H3.[cases$date>from & cases$date<to])
  }
  adjusted_reassortment = rbind(adjusted_reassortment, data.frame(virus = virus, date = date, method = method,
                                                                  peakCases = peakCases,
                                                                  removedSegments = removedSegments,
                                                                  cumulative = cumulative,
                                                                  median = median(rate), lower = HpdRate[1], upper = HpdRate[2]))
}



# plot cummulative vs. reassortment rates
subset = adjusted_reassortment
subset_rem = adjusted_reassortment
# subset_rem = data.frame()
# # subset_rem = subset[subset$removedSegments == "none",]
# # for each virus and year and method, keep only the removedSegments with the min
# # median rate
# for (v in unique(subset$virus)){
#   for (y in unique(format(subset$date, "%Y"))){
#     for (m in unique(subset$method)){
#       tmp = subset[subset$virus == v & format(subset$date, "%Y") == y & subset$method == m,]
#       if (nrow(tmp) > 1){
#         min_index = which.min(tmp$median)
#         subset_rem = rbind(subset_rem, tmp[min_index,])
#       }
#     }
#   }
# }

# for each season, use the year of $date -1 and add /$date$year
subset_rem$season = paste(as.numeric(format(subset_rem$date, "%y"))-1, "/", format(subset_rem$date, "%y"), sep = "")

# remove H1N1 season 2010/11
subset_rem = subset_rem[!(subset_rem$season == "21/22"),]
subset_rem = subset_rem[!(subset_rem$season == "22/23"),]

# calculate the correlation between the cumulative cases and the reassortment rate
# for each method and save it
corrdata=data.frame()
for (m in unique(subset_rem$method)){
  tmp = subset_rem[subset_rem$method == m,]
  corrdata = rbind(corrdata, data.frame(method = m, cor = cor(tmp$peakCases, tmp$median), p = cor.test(tmp$peakCases, tmp$median)$p.value))
}

# make the same plot, but with peak cases instead of cumulative cases
p1 = ggplot(subset_rem, aes(x = peakCases, y = median)) +
  geom_smooth(method = "lm", se = T, color="grey") +
  geom_point(aes(color=virus)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=virus)) +
  theme_minimal() +
  geom_text(data=corrdata, aes(label = paste("R^2=",round(cor,2), "\n p=", round(p,4)), x = 100, y = 0.4), hjust = 0, vjust = 1) +
  geom_text(data=subset_rem, aes(label = season,x = peakCases, y = median), hjust = 0, vjust = 1, size=2) +
  labs(x = "peak positive tests", y = "reassortment rate on extant lineages", title = "")+
  scale_color_manual(name="",values=c("H1N1"=h1n1color, "H3N2"=h3n2color))+
  facet_grid(method~.)+ theme(legend.position = "none")

# calculate the correlation between the cumulative cases and the reassortment rate
# for each method and save it
corrdata2=data.frame()
for (m in unique(subset_rem$method)){
  tmp = subset_rem[subset_rem$method == m,]
  corrdata2 = rbind(corrdata2, data.frame(method = m, cor = cor(tmp$cumulative, tmp$median), p = cor.test(tmp$cumulative, tmp$median)$p.value))
}

p2 = ggplot(subset_rem, aes(x = cumulative, y = median)) +
  geom_smooth(method = "lm", se = T, color="grey") +
  geom_point(aes(color=virus)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=virus)) +
  theme_minimal() +
  geom_text(data=corrdata2, aes(label = paste("R^2=",round(cor,2), "\n p=", round(p,4)), x = 100, y = 0.4), hjust = 0, vjust = 1) +
  geom_text(data=subset_rem, aes(label = season,x = cumulative, y = median), hjust = 0, vjust = 1, size=2) +
  labs(x = "cumulative positive tests", y = "reassortment rate on extant lineages", title = "")+
  scale_color_manual(name="",values=c("H1N1"=h1n1color, "H3N2"=h3n2color))+
  facet_grid(method~.) + theme(legend.position = "none")

legend = get_legend(
  ggplot(subset_rem, aes(x = peakCases, y = median, color = virus)) +
    geom_point() +
    scale_color_manual(name="", values=c("H1N1"=h1n1color, "H3N2"=h3n2color)) +
    theme(legend.position = "top")
)

# Combine plots with one legend on top
p3 = ggarrange(p1, p2, ncol = 2, nrow = 1)
p = ggarrange(legend, p3, ncol = 1, nrow = 2, heights = c(1, 5))
plot(p)

ggsave("./../../Figures/ExtantReassortmentRateVsCases.pdf", p, width = 8, height = 4)

# do the same, but read in the reassortment rates from the log files of all runs
# with constant in the name, use the column reassortmentRate.1
files <- list.files("./out/", pattern="*.log", full.names=TRUE)
# only keep the ones that say constant
files = files[grep("constant", files)]
reassortment_rates = data.frame()
for (i in seq(1, length(files))) {
  t <- read.table(files[i], header = TRUE, sep = "\t")
  # take a 10% burnin
  t = t[round(nrow(t)*0.8):nrow(t),]
  # get the virus and the year from the tree file
  virus = str_extract(basename(files[i]), "H1N1|H3N2")
  year = as.numeric(str_extract(basename(files[i]), "2015|2016|2017|2018|2019|2020|2021|2022|2023"))

  # convert the year to a date
  date = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  rate = t$reassortmentRate.1
  rate2 = t$reassortmentRate.2
  
  HpdRate = quantile(rate, c(0.025, 0.975))
  HpdRate2 = quantile(rate2, c(0.025, 0.975))
  
  # get the corresponding cummulative incidence for the same virus and year
  from = as.Date(paste(year, 6, 1, sep = "-"), format = "%Y-%m-%d")
  to = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  if (virus == "H1N1"){
    cumulative = max(cases$cumulativeH1N1[cases$date>from & cases$date<to])
    peakCases = max(who$cases[who$date>from & who$date<to & who$virus == "H1N1" & who$location=="North America"])
    peakCasesEU = max(who$cases[who$date>from & who$date<to & who$virus == "H1N1" & who$location=="Western Europe"])
  }else{
    cumulative = max(cases$cumulativeH3N2[cases$date>from & cases$date<to])
    peakCases = max(who$cases[who$date>from & who$date<to & who$virus == "H3N2" & who$location=="North America"])
    peakCasesEU = max(who$cases[who$date>from & who$date<to & who$virus == "H3N2" & who$location=="Western Europe"])
  }
  peakILI = max(ili$X..WEIGHTED.ILI[ili$date>from & ili$date<to])

  # get the logInfected.1,2,... with the largest median Ne
  logInfected = t[,grep("logInfected", colnames(t))][,2:6]
  maxNe = which.max(apply(logInfected, 2, median))
  Ne = exp(logInfected[,maxNe])
  reassortment_rates = rbind(reassortment_rates, data.frame(virus = virus, date = date, cumulative = cumulative, 
                                                              peakCases = peakCases,
                                                              peakCasesEU = peakCasesEU,
                                                              peakILI=peakILI,
                                                              median = median(rate), 
                                                              lower = HpdRate[1], 
                                                              upper = HpdRate[2],
                                                              median2 = median(rate2),
                                                              lower2 = HpdRate2[1],
                                                              upper2 = HpdRate2[2],
                                                              medianNe = median(Ne),
                                                              lowerNe = quantile(Ne, 0.025),
                                                              upperNe = quantile(Ne, 0.975)
                                                              ))
}


# make a new data frame that has arrows of length log2(peakILIratio$normalized)
# times some constant. The arrows should start at x=cumulative, y = median
# and then be horizontal
arrows = data.frame()
scale_for_arrow = 0.1
reassortment_rates$adjusted = NA
for (i in seq(1, length(reassortment_rates$date))){
  # get the corresponding peakILIratio
  peakILI = peakILIratio$normalized[peakILIratio$year == as.numeric(format(reassortment_rates$date[i], "%Y"))-1]
  arrows = rbind(arrows, data.frame(x = reassortment_rates$peakCases[i], y = reassortment_rates$median[i], 
                                     xend = reassortment_rates$peakCases[i] * peakILI, yend = reassortment_rates$median[i]))
  reassortment_rates$adjusted[i] = reassortment_rates$peakCases[i] * peakILI
}

reassortment_rates$year=format(reassortment_rates$date-365, "%Y")
# calculate the correlation between the cumulative cases and the reassortment rate
# for each method and save it
corrdata = data.frame(cor = cor(reassortment_rates$cumulative, reassortment_rates$median), 
                      p = cor.test(reassortment_rates$cumulative, reassortment_rates$median)$p.value)

# plot the reassortment rates from the log files
p = ggplot(reassortment_rates, aes(x = cumulative, y = median)) + 
  geom_segment(data=arrows, aes(x = x, y = y, xend = xend, yend = yend), alpha=0.2,arrow = arrow(length = unit(0.1, "inches")))+
  geom_point(aes(color=virus)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, color=virus)) +
  theme_minimal() + 
  geom_text(data=corrdata, aes(label = paste("R^2=",round(cor,2), "\n p=", round(p,6)), x = 100, y = 0.5), color="black", hjust = 0, vjust = 1) +
  geom_smooth(method = "lm", se=T, formula = y ~ x + 0) +
  labs(x = "cumulative CDC influenza incidence until May 31st", y = "reassortment rate", title = "")+
  geom_text(data=reassortment_rates, aes(label = format(date-365, "%Y"), x = cumulative, y = median), hjust = 0, vjust = 1) 
plot(p)


# make a new data frame with arrows for the below biases, make the
# arrows from peak cases to 1.2*peakcases at the median estimate
# season 2017/2018 higher H1N1 in europe but much lower H3N2
# season 2019/2020 much higher H3N2 in euope, similar H1N1
vals = reassortment_rates[(reassortment_rates$virus == "H3N2" & reassortment_rates$year == "2017"),]
bias_arrows = data.frame(x=vals$peakCases, y=vals$median, xend=vals$peakCases-500, yend=vals$median)
vals = reassortment_rates[(reassortment_rates$virus == "H1N1" & reassortment_rates$year == "2017"),]
bias_arrows = rbind(bias_arrows, 
                    data.frame(x=vals$peakCases, y=vals$median, xend=vals$peakCases+500, yend=vals$median))
vals = reassortment_rates[(reassortment_rates$virus == "H3N2" & reassortment_rates$year == "2019"),]
bias_arrows = rbind(bias_arrows, 
                    data.frame(x=vals$peakCases, y=vals$median, xend=vals$peakCases+500, yend=vals$median))



# reassortment_rates = reassortment_rates[!(reassortment_rates$virus == "H3N2" & reassortment_rates$year == "2019"),]
# reassortment_rates = reassortment_rates[!(reassortment_rates$virus == "H3N2" & reassortment_rates$year == "2016"),]
# reassortment_rates = reassortment_rates[!(reassortment_rates$year == "2022"),]
# reassortment_rates = reassortment_rates[!(reassortment_rates$year == "2021"),]

# make the same plot, but using adjusted instead of cumulative
corrdata = data.frame(cor = cor(reassortment_rates$peakCases, reassortment_rates$median), 
                      p = cor.test(reassortment_rates$peakCases, reassortment_rates$median)$p.value)
p = ggplot(reassortment_rates, aes(x = peakCases, y = median,color=virus)) + 
  geom_smooth(method = "lm", se=T, color="grey") +
  geom_point(aes()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  

  theme_minimal() + 
  geom_text(data=corrdata, aes(label = paste("R^2=",round(cor,2), "\n p=", round(p,6)), x = 100, y = 0.5), color="black", hjust = 0, vjust = 1) +
  labs(x = "peak cases during the season ", y = "reassortment rate", title = "")+
  scale_color_manual(name="",values=c("H1N1"=h1n1color, "H3N2"=h3n2color))+
  geom_text(data=reassortment_rates, aes(label = format(date-365, "%Y"), x = peakCases, y = median), hjust = 0, vjust = 1)
plot(p)

ggsave("./../../Figures/ReassortmentPeakCases.pdf", p, width = 6, height = 4)





p = ggplot(reassortment_rates, aes(x = median2, y = median,color=virus)) + 
  geom_smooth(method = "lm", se=T, color="grey") +
  geom_point(aes()) + 
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_errorbarh(aes(xmin=lower2, xmax=upper2)) +
  theme_minimal() + 
  # geom_text(data=corrdata, aes(label = paste("R^2=",round(cor,2), "\n p=", round(p,6)), x = 100, y = 0.5), color="black", hjust = 0, vjust = 1) +
  labs(x = "prev reassortment rates ", y = "reassortment rate", title = "")+
  scale_color_manual(name="",values=c("H1N1"=h1n1color, "H3N2"=h3n2color))+
  geom_text(data=reassortment_rates, aes(label = format(date-365, "%Y"), x = median2, y = median), hjust = 0, vjust = 1)
plot(p)




model <- glm(median ~ peakCases + peakCasesEU + 1, 
             data = reassortment_rates, 
             family = gaussian())
summary(model)
# Compute partial residuals for peakCases
residuals <- residuals(model)
partial_residuals <- residuals + coef(model)["peakCases"] * reassortment_rates$peakCases

# Create a data frame for plotting
plot_data <- data.frame(
  peakCases = reassortment_rates$peakCases,
  partial_residuals = partial_residuals
)

# Plot the partial residuals against peakCases
ggplot(plot_data, aes(x = peakCases, y = partial_residuals)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Partial Residual Plot for peakCases",
       x = "peakCases",
       y = "Partial Residuals (Corrected for peakCasesEU)") +
  theme_minimal()

corrdata = data.frame(cor = cor(reassortment_rates$peakCases, reassortment_rates$medianNe), 
                      p = cor.test(reassortment_rates$peakCases, reassortment_rates$medianNe)$p.value)
p = ggplot(reassortment_rates, aes(x = peakCases, y = medianNe)) + 
  geom_point(aes(color=virus)) + 
  geom_errorbar(aes(ymin=lowerNe, ymax=upperNe, color=virus)) +
  theme_minimal() + 
  geom_text(data=corrdata, aes(label = paste("R^2=",round(cor,2), "\n p=", round(p,4)), x = 100, y = 1000), color="black", hjust = 0, vjust = 1) +
  geom_smooth(method = "lm", se=T, formula = y ~ x + 0, color="red") +
  # geom_smooth(method = "lm", se=F, color="green") +
  # scale_y_log10() +
  # scale_x_log10() +
  coord_cartesian(ylim = c(0, 1000)) +
  labs(x = "CDC influenza seasonal incidence adjusted using ILI", y = "reassortment rate", title = "")+
  geom_text(data=reassortment_rates, aes(label = format(date-365, "%Y"), x = peakCases, y = medianNe), hjust = 0, vjust = 1)
plot(p)

ggsave("./../../Figures/NePeakCases.pdf", p, width = 8, height = 4)




# do the same, but read in the reassortment rates from the log files of all runs
# with constant in the name, use the column reassortmentRate.1
files <- list.files("./out/", pattern="*.log", full.names=TRUE)
# only keep the ones that say constant
files = files[grep("ne", files)]
reassortment_rates_ne = data.frame()
for (i in seq(1, length(files))) {
  t <- read.table(files[i], header = TRUE, sep = "\t")
  # take a 10% burnin
  t = t[round(nrow(t)*0.8):nrow(t),]
  # get the virus and the year from the tree file
  virus = str_extract(basename(files[i]), "H1N1|H3N2")
  year = as.numeric(str_extract(basename(files[i]), "2015|2016|2017|2018|2019|2020|2021|2022|2023"))
  
  # convert the year to a date
  date = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  rate = t$To
  HpdRate = quantile(rate, c(0.025, 0.975))
  
  # get the corresponding cummulative incidence for the same virus and year
  from = as.Date(paste(year, 6, 1, sep = "-"), format = "%Y-%m-%d")
  to = as.Date(paste(year+1, 5, 31, sep = "-"), format = "%Y-%m-%d")
  if (virus == "H1N1"){
    cumulative = max(cases$cumulativeH1N1[cases$date>from & cases$date<to])
    peakCases = max(who$cases[who$date>from & who$date<to & who$virus == "H1N1" & who$location=="North America"])
    peakCasesEU = max(who$cases[who$date>from & who$date<to & who$virus == "H1N1" & who$location=="Western Europe"])
  }else{
    cumulative = max(cases$cumulativeH3N2[cases$date>from & cases$date<to])
    peakCases = max(who$cases[who$date>from & who$date<to & who$virus == "H3N2" & who$location=="North America"])
    peakCasesEU = max(who$cases[who$date>from & who$date<to & who$virus == "H3N2" & who$location=="Western Europe"])
  }
  peakILI = max(ili$X..WEIGHTED.ILI[ili$date>from & ili$date<to])
  
  
  # get the logInfected.1,2,... with the largest median Ne
  logInfected = t[,grep("logInfected", colnames(t))][,2:6]
  maxNe = which.max(apply(logInfected, 2, median))
  Ne = exp(logInfected[,maxNe])
  reassortment_rates_ne = rbind(reassortment_rates_ne, data.frame(virus = virus, date = date, cumulative = cumulative, 
                                                            peakCases = peakCases,
                                                            peakCasesEU = peakCasesEU,
                                                            peakILI=peakILI,
                                                            median = median(rate), 
                                                            lower = HpdRate[1], 
                                                            upper = HpdRate[2],
                                                            medianNe = median(Ne),
                                                            lowerNe = quantile(Ne, 0.025),
                                                            upperNe = quantile(Ne, 0.975)
  ))
}



# make a function that returns the probability of observing a reassortment event
# given a binomial coefficient p the n0 segs left and no_segs right
binomial_coefficient = function(p, left, right){
  return(p^left*(1-p)^right + p^right*(1-p)^left)
}

# Function to compute the negative log-likelihood
neg_log_likelihood <- function(p, left_counts, right_counts) {
  log_likelihoods <- log(binomial_coefficient(p, left_counts, right_counts))
  return(-sum(log_likelihoods))
}

files <- list.files("./out/", pattern="*.events.reassortment", full.names=TRUE)
results = data.frame()
no_splits = data.frame()
for (file in files){
  t = read.table(file, header = FALSE, sep = "\t")
  # make two more columns that each count the elements in V2 or V3
  t$V3_count = sapply(t$V3, function(x) length(unlist(strsplit(x, ","))))
  t$V4_count = sapply(t$V4, function(x) length(unlist(strsplit(x, ","))))

  # split the absename on _ or dot
  fname = strsplit(basename(file), "[_.]")[[1]]
  for (iteration in unique(t$V1)){
    tmp = t[t$V1 == iteration & t$V2<0.8,]
    # get the mle estimate using the probs from binomial_coefficient for tmp

    # Get counts
    left_counts <- tmp$V3_count
    right_counts <- tmp$V4_count
    # Compute the MLE using optimization
    mle_result <- optim(par = 0.7, fn = neg_log_likelihood, left_counts = left_counts, right_counts = right_counts, method = "Brent", lower = 0.5, upper = 1)
    # Store the results
    results <- rbind(results, data.frame(virus = fname[1], year=fname[2], method=fname[3], iteration = iteration, mle = mle_result$par, count = length(tmp$V1)))
    
    for (i in seq(0, 6)){
      for (j in seq(i+1, 7)){
        # count how often there is i in V3 and j in V4 and the other way around
        overlap = sum(grepl(i, tmp$V3) & grepl(j, tmp$V4))
        overlap = overlap + sum(grepl(j, tmp$V3) & grepl(i, tmp$V4))
        
        # count the number of times they took the same path
        same = sum(grepl(i, tmp$V3) & grepl(j, tmp$V3))
        same = same + sum(grepl(i, tmp$V4) & grepl(j, tmp$V4))
        
        # save to no_splits
        no_splits = rbind(no_splits, data.frame(virus = fname[1], year=fname[2], method=fname[3], 
                                                iteration = iteration, i = i, j = j, 
                                                overlap = overlap, same = same))
        
        

      }
    }
    
  }
}


# plot the results
p = ggplot(results, aes(x = year, fill=virus, linetype=method, y=mle)) + 
  geom_violin()+
  theme_minimal() + 
  scale_fill_manual(values=c("H1N1"=h1n1color, "H3N2"=h3n2color))+
  # scale_color_ordinal()+
  labs(x = "MLE estimate of reassortment probability", y = "count", title = "")
plot(p)

# make a new data frame, and loop over all season and method, for each season, where there is both H1N1 and H3N2
# compute the median mle estimate over the posterior distribution for H1N1 and H3N2, then go to thenext
mle_data = data.frame()
for (y in unique(results$year)){
  for (m in unique(results$method)){
    tmp = results[results$year == y & results$method == m,]
    if (length(unique(tmp$virus)) == 2){
      mle_data = rbind(mle_data, data.frame(year = y, method = m, 
                                            h3n2 = median(tmp[tmp$virus == "H3N2",]$mle), 
                                            h1n1 = median(tmp[tmp$virus == "H1N1",]$mle),
                                            h3n2_lower = quantile(tmp[tmp$virus == "H3N2",]$mle, 0.025),
                                            h3n2_upper = quantile(tmp[tmp$virus == "H3N2",]$mle, 0.975),
                                            h1n1_lower = quantile(tmp[tmp$virus == "H1N1",]$mle, 0.025),
                                            h1n1_upper = quantile(tmp[tmp$virus == "H1N1",]$mle, 0.975)))
    }
  }
}

# plot the results and compute the correlation between the two
p = ggplot(mle_data, aes(x = h3n2, y = h1n1, color=method)) + 
  geom_errorbar(aes(ymin=h1n1_lower, ymax=h1n1_upper)) +
  geom_errorbarh(aes(xmin=h3n2_lower, xmax=h3n2_upper)) +
  geom_point() + 
  geom_smooth(method = "lm", se=T) +
  theme_minimal() + 
  labs(x = "H3N2", y = "H1N1", title = "")
plot(p)

# for each segment combination, plot the number of times they took the same path
p = ggplot(no_splits[no_splits$method=="constant", ], aes(x = interaction(virus, j, i), y = overlap/(same+overlap), fill=as.character(j))) + 
  geom_violin() + 
  theme_minimal() + 
  facet_grid(year~virus) +
  labs(x = "segment 1", y = "segment 2", fill = "overlap", title = "")
plot(p)

