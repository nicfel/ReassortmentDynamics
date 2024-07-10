library(stringr)
library(seqinr)
library(ggplot2)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

data = data.frame()
# get all the log files in /out/
log.files <- list.files("./out/", pattern="*.log", full.names=TRUE)
# only keep the files with constant in the name
log.files <- log.files[grep("constant", log.files)]
log.files <- log.files[-grep("wgs", log.files)]
# loop over all log files and read them in as a data frame
dynamics = data.frame()
timepoints = seq(0,5, length.out=15)
for (i in seq(1, length(log.files))) {
  # read in the log file
  t <- read.table(log.files[i], header=TRUE)
  # take a 10% burnin
  t <- t[-c(1:round(0.1*nrow(t))),]
  
  # get the name of the segments
  seg1 <- strsplit(log.files[i], "_")[[1]][2]
  tmp <- gsub(".log","", strsplit(log.files[i], "_")[[1]][[3]])
  seg2 = strsplit(tmp, "\\.")[[1]][1]
  method = strsplit(tmp, "\\.")[[1]][2]
  
  print(paste(seg1, seg2))
  # get the number of reassortment events, the reassortment rate and the median tree
  # length of both segments
  values = c("reassortmentRate", "network.reassortmentNodeCount", "seg1tree.treeLength", "seg2tree.treeLength")
  # for each, get the median and the 95% HPD and add them to a dataframe data as one row per log file
  data = rbind(data, data.frame(seg1=seg1, seg2=seg2, method=method, 
                                reassortmentRate=median(t$reassortmentRate.1),
                                reassortmentRate.l=quantile(t$reassortmentRate.1, 0.025),
                                reassortmentRate.u=quantile(t$reassortmentRate.1, 0.975),
                                network.reassortmentNodeCount=median(t$network.reassortmentNodeCount),
                                network.reassortmentNodeCount.l=quantile(t$network.reassortmentNodeCount, 0.025),
                                network.reassortmentNodeCount.u=quantile(t$network.reassortmentNodeCount, 0.975),
                                seg1tree.treeLength=median(t$seg1tree.treeLength),
                                seg1tree.treeLength.l=quantile(t$seg1tree.treeLength, 0.025),
                                seg1tree.treeLength.u=quantile(t$seg1tree.treeLength, 0.975),
                                seg2tree.treeLength=median(t$seg2tree.treeLength),
                                seg2tree.treeLength.l=quantile(t$seg2tree.treeLength, 0.025),
                                seg2tree.treeLength.u=quantile(t$seg2tree.treeLength, 0.975)
                                ))
  
  
  # If the name of the log file contains ne or variable, read in the log Ne
  if (grepl("ne", log.files[i]) | grepl("variable", log.files[i])) {
    # loop over every time point and get the Ne using column logInfected.timepoint
    for (timepoint in seq(1,length(timepoints))) {
      # get the column name
      colname = paste("logInfected.", timepoint, sep="")
      # if the log files is variable
      if (grepl("variable", log.files[i])) {
        colname2 = paste("InfectedToRho.", timepoint, sep="")
        vals = t[,colname2]
      }else{
        colname2 = "InfectedToRho"
        vals = log(t[,colname2])
      }

      # get the median and the 95% HPD of the Ne
      dynamics = rbind(dynamics, data.frame(seg1=seg1, seg2=seg2, 
                                            method=method, 
                                            timepoint=timepoints[timepoint],
                                            Ne=median(t[,colname]),
                                            Ne.l=quantile(t[,colname], 0.025),
                                            Ne.u=quantile(t[,colname], 0.975),
                                            prev=median(t[,colname]+vals),
                                            prev.l=quantile(t[,colname]+vals, 0.025),
                                            prev.u=quantile(t[,colname]+vals, 0.975)
                                            ))
      
    }
  }
}

# data is an upper triangular matrix, so we need to add the lower triangular part
data = rbind(data, data.frame(seg1=data$seg2, seg2=data$seg1, method=data$method, 
                              reassortmentRate=data$reassortmentRate,
                              reassortmentRate.l=data$reassortmentRate.l,
                              reassortmentRate.u=data$reassortmentRate.u,
                              network.reassortmentNodeCount=data$network.reassortmentNodeCount,
                              network.reassortmentNodeCount.l=data$network.reassortmentNodeCount.l,
                              network.reassortmentNodeCount.u=data$network.reassortmentNodeCount.u,
                              seg1tree.treeLength=data$seg2tree.treeLength,
                              seg1tree.treeLength.l=data$seg2tree.treeLength.l,
                              seg1tree.treeLength.u=data$seg2tree.treeLength.u,
                              seg2tree.treeLength=data$seg1tree.treeLength,
                              seg2tree.treeLength.l=data$seg1tree.treeLength.l,
                              seg2tree.treeLength.u=data$seg1tree.treeLength.u
                              ))

# plot the reassortment rates a matrix over seg1 and seg2
p = ggplot(data, aes(x=seg2, y=reassortmentRate,ymin = reassortmentRate.l, ymax = reassortmentRate.u, color=reassortmentRate)) + 
  geom_point()+
  geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Reassortment rate", x="Segment 1", y="Reassortment rate") +
  facet_grid(seg1~.) + scale_y_log10()
  theme_minimal()
plot(p)

# plot the means as a tile plot
p = ggplot(data, aes(x=seg2, y=seg1, fill=log(reassortmentRate))) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Reassortment rate", x="Segment 1", y="Segment 2") +
  scale_fill_viridis_c()+
  theme_minimal()
plot(p)

# plot the dynamics of Ne over time
# p = ggplot(dynamics[dynamics$timepoint<3,], aes(x=timepoint, y=Ne, ymin = Ne.l, ymax = Ne.u, color=seg1)) + 
#   geom_line()+
#   geom_ribbon(alpha=0.2, fill=NA)+
#   labs(title="Ne dynamics", x="Time", y="Ne") +
#   facet_grid(seg2~method) +
#   theme_minimal()
# plot(p)

# plot the dynamics of prevalence over time
# p = ggplot(dynamics[dynamics$timepoint<3 & dynamics$method=="variable",], 
#            aes(x=timepoint, y=exp(prev)/150, ymin = exp(prev.l)/150, ymax = exp(prev.u)/150, color=seg1, fill=seg1)) + 
#   geom_line()+
#   geom_ribbon(alpha=0.1)+
#   labs(title="Reassortment Rate", x="Time", y="Prevalence") +
#   facet_grid(seg1~seg2) +
#   coord_cartesian(ylim=c(0,0.01)) +
#   theme_minimal()
# plot(p)






