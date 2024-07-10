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

#define the number of states
nstates = 20

# read in the file SIR_simulations.txt that contains the simulated values
simulated = read.table("eventRateSIR_simulations.txt", header=TRUE, sep="\t")


print("read in SIR files")

# get the true values for the prevalence by readin in the files in /master
attackRate = data.frame()
# data frame to keep track of the timings of reassortment events to plot
reassortment.events = data.frame()
# get all .log files in .master
logfiles <- list.files("./master", pattern="*.log", full.names=TRUE)
# loop over all files
for (i in seq(1, length(logfiles))){
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
  # Do the same for lineages
  lins <- str_split(str_replace_all(data[[3]], "\\[|\\]", ""), ', ')[[1]]
  # get the total process time
  total_process_time = as.numeric(strsplit(sir[[length(sir)]], split=":")[[1]][[2]])

  events = 0
  
  # compute the number of lineages over time for each of nstates as a vector of 0's
  currlins = matrix(0, nrow=1, ncol=nstates)
  
  # get the state of the first lineage
  currlins = 1
  currtime = as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[2]]) # keeps the current time
  totalLength = 0;
  #also keep track of the mrca time
  mrca_time = -currtime
  
  for (j in seq(length(lins),1,-1)){
    # split the string on : and get the second group
    tmp = strsplit(lins[[j]], split=":")[[1]]
    # update the time since last reassortment
    totalLength = totalLength +
      currlins * (as.numeric(tmp[[2]]) - currtime)
    currtime = as.numeric(tmp[[2]])
    if (tmp[[1]] == "2"){ # coalescent event
      currlins = currlins + 1
    }else if (tmp[[1]] == "3"){ # migration event
    }else{ # sampling or co-infection event
      if (tmp[[1]] == "1"){
        events = events+1
      }
      currlins = currlins - 1
    }
    
    # check if any elements of currLins are below 0
    if (any(currlins<0)){
      print("negative lineages")
      error
    }
  }
  mrca_time = mrca_time + currtime

  # read in the xml and find the line with <populationSize spec="IntegerParameter" value="..."", save the populaiton sizes into a vector
  xmlfile = gsub(".log", ".xml", logfiles[i])
  xml = readLines(xmlfile)
  popSize = as.numeric(strsplit(strsplit(xml[grep("<populationSize", xml)], split="\"")[[1]][4], split=" ")[[1]])
  transmissionRate = as.numeric(strsplit(strsplit(xml[grep("<transmissionRate", xml)], split="\"")[[1]][4], split=" ")[[1]])[[1]]

  # Initialize I as a matrix with nstates rows and one column
  I <- matrix(0, nrow=1, ncol=nstates)
  # get the state of the first event
  tmp <- str_split(sir[1], ':')[[1]]
  I[as.numeric(strsplit(tmp[[4]], split="_")[[1]][[1]])+1] = 1
  time_I <- 0
  for (j in seq_along(sir)) {
      tmp <- str_split(sir[j], ':')[[1]]
      type = as.numeric(strsplit(tmp[[4]], split="_")[[1]][[1]])+1
      if (tmp[1] == '1') {
        I[type] = I[type]+1
      } else if (tmp[1] == '3') {
        to = as.numeric(strsplit(tmp[[4]], split="_")[[1]][[2]])+1
        I[to] = I[to]+1
        I[type] = I[type]-1
        # migration event changes the popoulation size, should probably be adapted....
        popSize[to] = popSize[to]+1
        popSize[type] = popSize[type]-1
      }
  }
  
  if (any(I/popSize>1)){
    dsa
  }
  # keep track of run and the prevalence
  print(mean(I/popSize))
  attackRate = rbind(attackRate, data.frame(run=run, transmissionRate=transmissionRate, events=events, length=totalLength, 
                                            mrca_time=mrca_time,
                                            total_process_time=total_process_time,
                                            attackRate=mean(I/popSize)))
  
}

# plot the attack rate vs. the number of reassortment events
p = ggplot(attackRate, aes(y=events/length/transmissionRate, x=attackRate/mrca_time)) + 
  geom_point() + geom_smooth(method="lm") +
  # add a diagonal line
  geom_abline(intercept = 0, slope = 1, color="red") +
  # scale_y_log10() +
  # scale_x_log10() +
  theme_minimal() + 
  ylab("Reassortment rate") + 
  xlab("Attack rate") + 
  ggtitle("Attack rate vs. number of reassortment events") + 
  theme(plot.title = element_text(hjust = 0.5))
plot(p)

p = ggplot(attackRate, aes(x=events/transmissionRate, y=attackRate/mrca_time)) + 
  geom_point() + geom_smooth(method="lm") +
  # scale_y_log10() +
  # scale_x_log10() +
  theme_minimal() + 
  ggtitle("Attack rate vs. number of reassortment events") + 
  theme(plot.title = element_text(hjust = 0.5))
plot(p)



