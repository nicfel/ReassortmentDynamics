library(stringr)
library(seqinr)
library(ggplot2)
library("colorblindr")
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

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

segments = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)



# get all tree files with rep0 in out
tree.files = list.files("out", pattern="rep0.trees", full.names=T)
# combine the tree file with the rep1 and rep2 using log combiner
for (i in 1:length(tree.files)) {
  # system(
  #        paste("/Applications/BEAST\\ 2.7.6/bin/logcombiner -burnin 80 -log", tree.files[i],
  #              gsub("rep0", "rep1", tree.files[i]), gsub("rep0", "rep2", tree.files[i]),
  #              "-o",
  #              gsub(".rep0.trees", ".combined.trees", tree.files[i])))
}
# get all the combined trees
tree.files = list.files("out", pattern="combined.trees", full.names=T)
# run them through the reassortmentRtaeOverTime app
for (i in 1:length(tree.files)) {
  # system(
  #        paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentNetworkSummarizer",
  #              "-burnin 0 -positions none",
  #              "-followSegment 0",tree.files[i],
  #              gsub(".trees", ".tree", tree.files[i])))
  for (a in seq(0,7)){
    for (b in seq(a+1,7)){
      remove_segments = seq(0,7)
      # remove the entries ==a and ==b
      remove_segments = remove_segments[remove_segments!=a]
      remove_segments = remove_segments[remove_segments!=b]
      # system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
      #        paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentNetworkSummarizer",
      #                        "-burnin 0 -positions none -removeSegments",paste(remove_segments, collapse=","),
      #                        "-followSegment", a, tree.files[i],
      #                        gsub(".trees", paste(".", a,"_",b, ".tree", sep=""), tree.files[i])))
      # 
      # system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
      #        paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentOverTime",
      #              "-burnin 0 -removeSegments",paste(remove_segments, collapse=",") , tree.files[i],
      #              gsub(".trees", paste(".", a,"_",b, ".reassortment.txt", sep=""), tree.files[i])))
      
    }
  }
}






data = data.frame()

burnin = 0.8

# get all files ending in *varying.log files in out
t1 = read.table("out/H5N1_wgs_proportional.varying.rep0.log", header = TRUE, sep = "\t")
t2 = read.table("out/H5N1_wgs_proportional.varying.rep1.log", header = TRUE, sep = "\t")
t3 = read.table("out/H5N1_wgs_proportional.varying.rep1.log", header = TRUE, sep = "\t")
# combined after 10% burnin
t = rbind(t1[round(nrow(t1)*burnin):nrow(t1),], t2[round(nrow(t2)*burnin):nrow(t2),], t3[round(nrow(t3)*burnin):nrow(t3),])
# read in the xml file xmls/H5N1_wgs.varying.rep0.xml
xmlFile = readLines("xmls/H5N1_wgs_proportional.varying.rep0.xml")
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
# clear xmlFIle from memory
rm(xmlFile)
# make two datframes, the first reassortment, reads in the piecewise constant reassortment rates using rateShifts2
reassortment = data.frame()
timepoints = as.numeric(strsplit(rateShifts2, split=" ")[[1]])
for (i in 1:(length(timepoints))) {
  rate = t[, paste("reassortmentRate", i, sep=".")]
  reassortment = rbind(reassortment, data.frame(from=mrsi-timepoints[i]*365, 
                                                rate=median(rate), 
                                                upper.5=quantile(rate, 0.75), lower.5=quantile(rate, 0.25),
                                                upper=quantile(rate, 0.975), lower=quantile(rate, 0.025)))
  reassortment = rbind(reassortment, data.frame(from=mrsi-timepoints[i+1]*364, 
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
p = ggplot(ne, aes(x=time)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) + 
  geom_ribbon(aes(ymin=lower.5, ymax=upper.5), alpha=0.75) + 
  # fill in the area between lower and upper for the reassortment rates
  geom_ribbon(data=reassortment, aes(x=from, ymin=lower.5*100, ymax=upper.5*100), alpha=0.85, fill="#ff7f00") +
  geom_ribbon(data=reassortment, aes(x=from, ymin=lower*100, ymax=upper*100), alpha=0.5, fill="#ff7f00") +
  # add a second axis for the reassortment rates
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "Reassortment rate")) +
  coord_cartesian(xlim=c(mrsi-365*2.5, mrsi), ylim=c(0,250)) +
  xlab("") + ylab("Effective population size")+
  theme_minimal()
plot(p)
ggsave("./../../Figures/H5N1_constant_dynamics.pdf", p, width = 9, height = 4)


# read in the constant files
t1 = read.table("out/H5N1_wgs_proportional.variable.rep0.log", header = TRUE, sep = "\t")
t2 = read.table("out/H5N1_wgs_proportional.variable.rep1.log", header = TRUE, sep = "\t")
t3 = read.table("out/H5N1_wgs_proportional.variable.rep2.log", header = TRUE, sep = "\t")
# combined after 10% burnin
t = rbind(t1[round(nrow(t1)*burnin):nrow(t1),], t2[round(nrow(t2)*burnin):nrow(t2),], t3[round(nrow(t3)*burnin):nrow(t3),])
# read in the xml file xmls/H5N1_wgs.varying.rep0.xml
xmlFile = readLines("xmls/H5N1_wgs_proportional.constant.rep0.xml")
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
mdsi = min(as.Date(dates))-20/12*365

# get the line with <rateShifts id="rateShifts2" spec="RealParameter" value="
rateShifts2 = strsplit(xmlFile[grep("<rateShifts id=\"rateShifts2\" spec=\"RealParameter\" value=\"", xmlFile)], 
                       split="\"")[[1]][6]
# clear xmlFIle from memory
rm(xmlFile)

maxrate = max(as.numeric(strsplit(rateShifts, " ")[[1]]))
timestep=maxrate/1000
# make two datframes, the first reassortment, reads in the piecewise constant reassortment rates using rateShifts2
reassortment_variable = data.frame()
currtime = mrsi
currint = 0
while (currtime >mdsi) {
  rate = t[, paste("reassortment", currint, sep="")]
  reassortment_variable = rbind(reassortment_variable, data.frame(from=currtime, 
                                                upper.5=quantile(rate, 0.75), 
                                                lower.5=quantile(rate, 0.25),
                                                upper=quantile(rate, 0.975), 
                                                lower=quantile(rate, 0.025)))
  currtime = mrsi - timestep*currint*365
  currint = currint + 1
  
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
# # if from is NA in reassortment_constant, put the min time of ne_constant$time there
reassortment_variable$from[is.na(reassortment_variable$from)] = min(ne$time)


# plot the two
p = ggplot(ne, aes(x=time)) + 
  geom_ribbon(aes(ymin=lower.5, ymax=upper.5), alpha=1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
  geom_ribbon(data=reassortment_variable, aes(x=from, ymin=lower*100, ymax=upper*100), alpha=0.5, fill="#ff7f00") +
  geom_ribbon(data=reassortment_variable, aes(x=from, ymin=lower.5*100, ymax=upper.5*100), alpha=0.85, fill="#ff7f00") +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "Reassortment rate")) +
  coord_cartesian(xlim=c(mrsi-365*2.5, mrsi), ylim=c(0,250)) +
  xlab("") +ylab("Effective population size")+
  theme_minimal()
plot(p)
ggsave("./../../Figures/H5N1_variable_dynamics.pdf", p, width = 9, height = 4)



# read in the reassortment.txt files, the first row is the run number, the second
# is the time, the third is 1/rate
reassortment = data.frame()
maxtime = 3
for (i in 1:length(tree.files)) {
  t = read.table(gsub(".trees", ".reassortment.txt", tree.files[i]), header=F, sep="\t")
  for (run in unique(t$V1)) {
    tmp = t[t$V1==run,]
    timeFrom = c(0, tmp$V2)
    name = strsplit(tree.files[i], split="\\.")[[1]]
    reassortment = rbind(reassortment, data.frame(run=run,
                                                  from = timeFrom[1:(length(timeFrom)-1)],
                                                  to = tmp$V2,
                                                  rate=1/tmp$V3, 
                                                  file=name[[2]]))
  }
}

reassortmentfiles = list.files("out", pattern="reassortment.txt", full.names=T)
# remove the files with 7_7 and 7_8
reassortmentfiles = reassortmentfiles[!grepl("7_7", reassortmentfiles)]
reassortmentfiles = reassortmentfiles[!grepl("7_8", reassortmentfiles)]
counts = data.frame()
for (i in 1:length(reassortmentfiles)) {
  name = strsplit(reassortmentfiles[i], split="\\.")[[1]]
  if (length(name)==6){
    t = read.table(reassortmentfiles[i], header=F, sep="\t")
    fromto = strsplit(name[[4]], split="_")[[1]]
    # for each unique V0, compute how many instances of V2 <maxtime there are
    for (run in unique(t$V1)) {
      tmp = t[t$V1==run,]
      counts = rbind(counts, data.frame(run=run, 
                                       from=fromto[1], to=fromto[2], 
                                       count=sum(tmp$V2 < maxtime),
                                       file=name[[2]]))
      counts = rbind(counts, data.frame(run=run, 
                                        from=fromto[2], to=fromto[1], 
                                        count=sum(tmp$V2 < maxtime),
                                        file=name[[2]]))
      
    }
  }
}

segments = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

colors = c('#fdbf6f', '#e31a1c', '#ff7f00', '#a6cee3', '#1f78b4', '#fb9a99', '#b2df8a', '#33a02c')

colors = c('HA'='#fdbf6f', 'NA'='#ff7f00','MP'='#a6cee3', 'NP'='#1f78b4', 'PA'='#66c2a4','PB1'='#2ca25f','PB2'='#006d2c', 'NS'='#fb9a99')
# plot a histogram for file=="constant" for the count for each combination of from to
#rename from and to from "0" "1".. "7" to segments
counts$froms = segments[as.numeric(counts$from)+1]
counts$tos = segments[as.numeric(counts$to)+1]
counts$froms = factor(counts$froms, levels=labels(colors))
counts$tos = factor(counts$tos, levels=labels(colors))

p = ggplot(counts, aes(x=count, fill=tos, color=tos)) + 
  geom_density(alpha=0.4) + 
  xlab("Reassortment events") + 
  # remove y axis ticks anf labels
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  facet_grid(froms~file ) +
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  # ggtitle("Reassortment counts") + 
  theme_minimal()
plot(p)
ggsave("./../../Figures/H5N1_CoReassortmentDistr.pdf", p, width = 9, height = 4)


p <- ggplot(counts[counts$froms == "PB2" & counts$file == "varying",], aes(x = count, fill = tos, color = tos)) + 
  geom_density(alpha = 0.4) + 
  xlab("Reassortment events") + 
  scale_color_manual(name="", values = colors) +
  scale_fill_manual(name="", values = colors) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 40)) +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(), 
    axis.line.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top"
  )
plot(p)
ggsave("./../../Figures/H5N1_PB2_varying_CoReassortmentDistr.pdf", p, width = 6, height = 3)

# plot the mean counts as a tile plot
library(dplyr)
mean_counts <- counts %>%
  group_by(from, to, file) %>%
  summarise(mean_count = mean(count, na.rm = TRUE)) %>%
  ungroup()

mean_counts$froms = segments[as.numeric(mean_counts$from)+1]
mean_counts$tos = segments[as.numeric(mean_counts$to)+1]
mean_counts$froms = factor(mean_counts$froms, levels=labels(colors))
mean_counts$tos = factor(mean_counts$tos, levels=rev(labels(colors)))

p = ggplot(mean_counts, aes(x=froms, y=tos, fill=log(mean_count))) + 
  geom_tile() + facet_wrap(~file) + 
  scale_fill_viridis_c(name="log number\nevents")+
  xlab("") + 
  ylab("") + 
  ggtitle("") + theme_minimal()
plot(p)
ggsave("./../../Figures/H5N1_CoReassortmentMean.pdf", p, width = 9, height = 4)

# for each file in reassormtne compute a moving average from time =0 to time = 3
timepoints = seq(0, 7, 0.1)
rea_avg=data.frame()
for (f in unique(reassortment$file)) {
  tmp = reassortment[reassortment$file==f,]
  for (time in timepoints) {
    rate = tmp[tmp$from <= time & tmp$to >= time, "rate"]
    rea_avg = rbind(rea_avg, data.frame(from=time, to=time+0.1, rate=mean(rate), 
                                        lower=quantile(rate, 0.025), upper=quantile(rate, 0.975),
                                        file=f))
  }
}

# plot everything as a step function in ggplot
p = ggplot(rea_avg, aes(x=from, xend=to, y=rate, yend=rate)) + 
  geom_segment() +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
  facet_wrap(~file) +
  xlab("Time (years)") + 
  coord_cartesian(ylim=c(0, 2)) +
  ylab("Reassortment rate") + ggtitle("Reassortment rate over time") + theme_minimal()
plot(p)

