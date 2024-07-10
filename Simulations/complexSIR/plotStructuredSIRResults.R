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
const.est.data = data.frame()
# data frame to save the reassortment events
reassortment.data = data.frame()

#define the number of states
nstates = 20

# read in the file SIR_simulations.txt that contains the simulated values
simulated = read.table("structuredSIR_simulations.txt", header=TRUE, sep="\t")

# get all .log files in /out
logfiles <- list.files("./out", pattern="*.log", full.names=TRUE)
# remove all files that do not contain the word variable not variable
# logfiles = logfiles[grepl("variable", logfiles) ]

# remove all log files that are not run 13
# logfiles = logfiles[grepl("_13.", logfiles) ]

# loop over all files
for (i in seq(1, length(logfiles))){
# for (i in seq(1, 30)){
  # read in the log file to get the network height
  filename = gsub("./out/", "./master/", logfiles[[i]])  
  filename = gsub(".constant.", ".",filename)  
  filename = gsub(".ne.", ".", filename)  
  filename = gsub(".variable.", ".", filename)  
  print(filename)
  
  lines = readLines(filename)
  data = strsplit(lines[[2]], split="\t")[[1]]
  lins <- str_split(str_replace_all(data[[3]], "\\[|\\]", ""), ', ')[[1]]
  
  timediff = as.numeric(strsplit(lins[[1]], split=":")[[1]][[2]]) - 
    as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[2]])

  # set the time points for the parameters
  time_points = seq(0,timediff*1.1, length.out=15)
  # get 500 points between 0 and 16
  time_points2 = seq(0, timediff*1.1, length.out=500)
  
  t = read.table(logfiles[i], header=TRUE, sep="\t")
  # remove 10 % as burnin
  t = t[round(0.2*nrow(t)):nrow(t),]
  # calculate the ESS for columns 2 to 5
  ess = sapply(t[,2:5], function(x) effectiveSize(x))
  # calculate the minimum ESS
  minESS = min(ess)
  # if the minimum ESS is below 100, print the filename
  if (minESS < 0) {
    print(logfiles[i])
    next
  }
  # split on . and get the second group as the method used
  method = strsplit(logfiles[i], "\\.")[[1]][3]
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][4])
  # if the method is not constant
  if (!grepl("constant", method)){
    # initialize the prevalence vector
    prevalence = c()
    prevalencel = c()
    prevalenceu = c()
    time = c()
    # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
    for (j in seq(0, length(time_points2)-1)){
      # get the prevalence at the time point, if the label exists
      if (paste0("reassortment", j) %in% colnames(t)){
        prevalence = c(prevalence, median(t[,paste0("reassortment", j)]))
        # get the 2.5 and 97.5 quantile
        prevalencel = c(prevalencel, quantile(t[,paste0("reassortment", j)], 0.025))
        prevalenceu = c(prevalenceu, quantile(t[,paste0("reassortment", j)], 0.975))
        # time
        time = c(time, time_points2[j+1])
      }    
    }
    # get the value in simulated in the run row and "transmission" column
    transmision0 = simulated[run, "transmission"]
    # devide prevalence by transmission0
    prevalence = prevalence/transmision0
    prevalencel = prevalencel/transmision0
    prevalenceu = prevalenceu/transmision0
    # add the data to a data frame
    est.data = rbind(est.data, data.frame(method=method, run=run, time=time, 
            prevalence=prevalence, prevalencel=prevalencel, prevalenceu=prevalenceu)) 
  }else{
    # if the method is constant
    # get the value in simulated in the run row and "transmission" column
    transmision0 = simulated[run, "transmission"]
    # get the reassortment Rate
    prevalence = median(t[, "reassortmentRate"])
    # get the 2.5 and 97.5 quantile
    prevalencel = quantile(t[, "reassortmentRate"], 0.025)
    prevalenceu = quantile(t[, "reassortmentRate"], 0.975)
    # dubblicate the above vectors to have the same value two times
    prevalence = c(prevalence, prevalence)
    prevalencel = c(prevalencel, prevalencel)
    prevalenceu = c(prevalenceu, prevalenceu)
    # time
    # add to the data.frame
    const.est.data = rbind(const.est.data, data.frame(method=method, run=run, time=NA, 
            prevalence=prevalence/transmision0, prevalencel=prevalencel/transmision0, 
            prevalenceu=prevalenceu/transmision0))
  }
  reassortment.data = rbind(reassortment.data, data.frame(method=method, run=run, 
          events=median(t[, "network.reassortmentNodeCount"]), 
          eventsl=quantile(t[, "network.reassortmentNodeCount"], 0.025), 
          eventsu=quantile(t[, "network.reassortmentNodeCount"], 0.975)))
}
est.data$mrsi = NA
reassortment.data$true = NA
print("read in SIR files")

const.est.data$weightedPrevalence = NA

# get the true values for the prevalence by readin in the files in /master
true.data = data.frame()
# data frame to keep track of the instantaneous reassortment rate
inst.rate = data.frame()
# data frame to keep track of the timings of reassortment events to plot
reassortment.events = data.frame()
# get all .log files in .master
logfiles <- list.files("./master", pattern="*.log", full.names=TRUE)
# remove all log files that are not run 13
# logfiles = logfiles[grepl("_13.", logfiles) ]

# loop over all files
for (i in seq(1, length(logfiles))){
# for (i in seq(1, 10)){
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][4])
  
  # check if the run is in reassortment.data, if not, continue
  if (!run %in% reassortment.data$run){
    next
  }
  
  
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
  
  # compute the number of lineages over time for each of nstates as a vector of 0's
  currlins = matrix(0, nrow=1, ncol=nstates)
  
  # get the state of the first lineage
  currloc = as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[4]])
  currlins[currloc+1] = 1
  
  # get the state of the last lineage
  currtime = as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[2]]) # keeps the current time
  roottime = as.numeric(strsplit(lins[[1]], split=":")[[1]][[2]])
  lineages_frame = data.frame() # dataframe for calcs after
  # keep track of the product of the number of lineages * delta t
  weighted_time_since_reassortment = 0
  for (j in seq(length(lins),1,-1)){
    # split the string on : and get the second group
    tmp = strsplit(lins[[j]], split=":")[[1]]
    # add to the lineages data frame
    lineages_frame = rbind(lineages_frame, data.frame(start=currtime, end=as.numeric(tmp[[2]]), l=currlins))
    # update the time since last reassortment
    weighted_time_since_reassortment = weighted_time_since_reassortment +
      sum(currlins) * (as.numeric(tmp[[2]]) - currtime)
    
    
    currtime = as.numeric(tmp[[2]])
    if (tmp[[1]] == "2"){ # coalescent event
      currloc = as.numeric(tmp[[4]])+1
      currlins[currloc] = currlins[currloc] + 1
    }else if (tmp[[1]] == "3"){ # migration event
      from = as.numeric(strsplit(tmp[[4]], split="_")[[1]][[2]])+1
      to = as.numeric(strsplit(tmp[[4]], split="_")[[1]][[1]])+1
      currlins[from] = currlins[from] - 1
      currlins[to] = currlins[to] + 1
    }else{ # sampling or co-infection event
      if (tmp[[1]] == "1"){
        reassortment.events = rbind(reassortment.events, data.frame(run=run, time=as.numeric(tmp[[2]])))
        inst.rate = rbind(inst.rate, data.frame(run=run, time=as.numeric(tmp[[2]]), rate=weighted_time_since_reassortment))
        weighted_time_since_reassortment = 0;
      }
      currloc = as.numeric(tmp[[4]])+1
      currlins[currloc] = currlins[currloc] - 1
    }
    
    # check if any elements of currLins are below 0
    if (any(currlins<0)){
      print("negative lineages")
      error
    }
  }
  # 
  mrsi = strsplit(lins[[1]], split=":")[[1]][[2]]
  # get the final time
  final = strsplit(sir[[length(sir)]], split=":")[[1]][[2]]
  # make 200 points between 0 and final
  time_points = seq(0, final, length.out=200)
  # keeps track of which time point we are at
  tc = 2;

  # add the mrsi to the est.data
  est.data[est.data$run==run, "mrsi"] = as.numeric(mrsi)
  
  # Initialize I as a matrix with nstates rows and one column
  I <- matrix(0, nrow=1, ncol=nstates)

  # get the state of the first event
  tmp <- str_split(sir[1], ':')[[1]]
  I[as.numeric(strsplit(tmp[[4]], split="_")[[1]][[1]])+1] = 1
  # Icurrent keeps track of the current state
  Icurrent = I

  time_I <- 0
  # Keep also track of the co-infection event times
  co_inf_times <- numeric()  
  
  # read in the xml and find the line with <populationSize spec="IntegerParameter" value="..."", save the populaiton sizes into a vector
  xmlfile = gsub(".log", ".xml", logfiles[i])
  xml = readLines(xmlfile)
  popSizeCurrent = as.numeric(strsplit(strsplit(xml[grep("<populationSize", xml)], split="\"")[[1]][4], split=" ")[[1]])
  popSize = popSizeCurrent

  for (j in seq_along(sir)) {
      tmp <- str_split(sir[j], ':')[[1]]
      type = as.numeric(strsplit(tmp[[4]], split="_")[[1]][[1]])+1
      if (tmp[1] == '0') {
        # copy the last row and add it to the matrix
      } else if (tmp[1] == '1') {
        Icurrent[type] = Icurrent[type]+1
      } else if (tmp[1] == '2') {
        Icurrent[type] = Icurrent[type]-1
      } else if (tmp[1] == '3') {
        to = as.numeric(strsplit(tmp[[4]], split="_")[[1]][[2]])+1
        Icurrent[type] = Icurrent[type]-1
        Icurrent[to] = Icurrent[to]+1
        # migration event changes the popoulation size, should probably be adapted....
        popSizeCurrent[type] = popSizeCurrent[type]-1
        popSizeCurrent[to] = popSizeCurrent[to]+1
      }

      # save every 20 state to I and time_I
      # if (as.numeric(tmp[2])>=time_points[[tc]]){
        I = rbind(I, Icurrent)
        popSize = rbind(popSize, popSizeCurrent)
        time_I = c(time_I, as.numeric(tmp[2]))
      #   tc=tc+1
      # }        
  }


  
  # for each start and end time in lineages_frame, compute the prevalence, 
  # waited by the number of lineages in a given location and the time.
  # keep track of the total time*lineages to normalize in the end
  normalization_factor = 0
  tot_prevalence = 0
  for (j in seq(2, length(lineages_frame$start))){
    # get all elements in time_I that are between start and end to get the true number of infected for that time period
    tmp = time_I[time_I >= lineages_frame$start[j] & time_I <= lineages_frame$end[j]]
    # multiply the number of infected individuals in each state by the number of lineages
    # in the same state for lineages_frame[j,3:end], but in R language
    infected_state = I[time_I >= lineages_frame$start[j] & time_I < lineages_frame$end[j],]
    # multiply the matrix infected_state with the dataframe vector lins_state,
    # but fix the orientation of lins_state
    popSizeState = popSize[time_I >= lineages_frame$start[j] & time_I < lineages_frame$end[j],]
    ratio = infected_state / popSizeState
    
    lins_state <- as.matrix(lineages_frame[j, 3:ncol(lineages_frame)])
    # divide the number of infected individuals in each state by the population size
    # in the same state
    weighted_prevalence = c()
    if (length(ratio)>nstates){
      for (row in seq(1, length(ratio)/nstates)){
        weighted_prevalence[row] = sum(ratio[row,] * lins_state)
      }
      weighted_prevalence = weighted_prevalence/sum(lins_state)
    }else if (length(ratio)==nstates){
      weighted_prevalence = sum(ratio * lins_state)
      weighted_prevalence = weighted_prevalence/sum(lins_state)
    }


    # check if any elemente in weighted_prevalence is below 0
    if (any(weighted_prevalence<0)){
      print("negative prevalence")
      dsa
    }

    # now also weigh for the time of each interval
    total_time = tmp[length(tmp)] - tmp[1]
    # weight the prevalence_times_lineages by the time of each interval
    if (total_time>0){
      time_weighted_prevalence = sum(weighted_prevalence * diff(tmp))/total_time
      tot_prevalence = tot_prevalence + time_weighted_prevalence * total_time * sum(lins_state)
      # compute the normalization factor to later normalize by the total time
      normalization_factor = normalization_factor + sum(lins_state) * sum(diff(tmp))
    }else{
      time_weighted_prevalence = mean(weighted_prevalence);
    }
    
    if (length(tmp)>1){
      true.data = rbind(true.data, data.frame(run=run, time=(tmp[length(tmp)]+tmp[1])/2, 
                                              weighted_prev=time_weighted_prevalence,
                                              mean_prev=mean(ratio),
                                              max_prev= max(ratio),
                                              num_lineages=sum(lins_state),
                                              method="time varying"))
    }
    else{
      true.data = rbind(true.data, data.frame(run=run, time=tmp[1], 
                                              weighted_prev=time_weighted_prevalence,
                                              mean_prev=mean(ratio),
                                              max_prev= max(ratio),
                                              num_lineages=sum(lins_state),
                                              method="time varying"))
    }
  }
  
  print(tot_prevalence/normalization_factor)
  # add the weighted prevalence to the const.est.data$weightedPrevalence by matching runs
  const.est.data[const.est.data$run==run, "weightedPrevalence"] = tot_prevalence/normalization_factor

  # experienced prevalence for constant calculation experienced_prev=tot_prevalence/normalization_factor)
  # get the true number of reassortment events in this run by counting
  # how often the strsplit(lins, split=":")[[:]][[1]] starts with 1:
  reassortment.data[reassortment.data$run==run, "true"] = sum(sapply(strsplit(lins, split=":"), function(x) x[[1]])=="1")

  # # add the data to the data frame state by state
  # for (j in seq(1, nstates)){
  #   true.data = rbind(true.data, data.frame(run=run, time=time_I, infected=I[,j], sumPop=sum(popSize), popSize=popSize[[j]], state=j, experienced_prev=tot_prevalence/normalization_factor))
  # }
}


# make a manual color scale, using alternating black and white for factor(state) and red for factor(method)
# get the number of states
nstates = length(unique(true.data$state))
# get the number of methods
nmethods = length(unique(est.data$method))
colors = c(rep(c("grey", "#343232"), nstates/2),
            rep(c("red", "green"), nmethods))

# compute the population size weighted average prevalence for each method and run at 
# every time point across all the states for true.data
# avg.data = data.frame()
# for (i in seq(1, length(unique(true.data$run)))){
#   print(i)
#   tp = unique(true.data[true.data$run==unique(true.data$run)[[i]], "time"])
#   for (j in seq(1, length(tp))){
#     indices = true.data$run==unique(true.data$run)[[i]] & true.data$time==tp[j]
#     m.val.pS = mean(true.data[indices, "infected"] / true.data[indices, "popSize"])
#     m.val.i = sum(true.data[indices, "infected"] / true.data[indices, "popSize"] * true.data[indices, "infected"])/sum(true.data[indices, "infected"])
#     avg.data = rbind(avg.data, data.frame(run=unique(true.data$run)[[i]], time=tp[j], prevalence=m.val.pS, prevalence.i=m.val.i))
#   }
# }
reassortment.events$y = NA
for (run in unique(reassortment.events$run)){
  # find the maximum upper value of prevalenceu for the corresponding run
  maxu = max(est.data[est.data$run==run & (est.data$method=="variable"), "prevalenceu"])
  reassortment.events[reassortment.events$run==run, "y"] = maxu
}


runnr = unique(reassortment.data[reassortment.data$true>30 & !is.na(reassortment.data$true), "run"])

# only pick runnr for which true.data has method=="ne" for this run
use.runs = c()
for (run in runnr){
  if (sum(est.data$run==run & est.data$method=="variable" )>0){
    use.runs = c(use.runs, run)
  }
}


use.runs = use.runs[seq(1, 10)]
# plot the results, for true.data, plot the stacked prevalence over time
p = ggplot(true.data[is.element(true.data$run,use.runs),])+
      geom_ribbon(data=est.data[is.element(est.data$run,use.runs) & est.data$method=="variable",], aes(x=mrsi-time, ymin=prevalencel, ymax=prevalenceu, fill=method), alpha=0.25) +
      geom_line(data=est.data[is.element(est.data$run,use.runs) & est.data$method=="variable",], aes(x=mrsi-time, y=prevalence, color=method)) +
      geom_ribbon(data=est.data[is.element(est.data$run,use.runs) & est.data$method=="ne",], aes(x=mrsi-time, ymin=prevalencel, ymax=prevalenceu, fill=method), alpha=0.25) +
      geom_line(data=est.data[is.element(est.data$run,use.runs) & est.data$method=="ne",], aes(x=mrsi-time, y=prevalence, color=method)) +
      geom_line(aes(x=time, y=weighted_prev, color="lineage weighted mean prevalence across states"))+
      geom_line(aes(x=time, y=mean_prev, color="mean prevalence across states"))+
      # geom_line(aes(x=time, y=num_lineages/500, color="lineages through time"))+
      # geom_line(aes(x=time, y=max_prev, color="max prevalence across states"))+
      geom_point(data=reassortment.events[is.element(reassortment.events$run,use.runs),], aes(x=time, y=0.19), alpha=0.2, size=2)+
      # geom_ribbon(data=est.data[is.element(est.data$run,use.runs) & est.data$method=="ne",], aes(x=mrsi-time, ymin=prevalencel, ymax=prevalenceu, fill=method), alpha=0.25) +
      # geom_line(data=est.data[is.element(est.data$run,use.runs) & est.data$method=="ne",], aes(x=mrsi-time, y=prevalence, color=method)) +
      # plot the reassortment events as vertical lines
      # geom_density(data = inst.rate[is.element(inst.rate$run,use.runs),], aes(x=time, color="instantaneous reassortment rate"))+
      facet_wrap(~run, ncol=3, scales = "free") +
      theme_minimal() +
      coord_cartesian(ylim=c(0,0.2)) +
      scale_color_manual(values=c("#946853", "#000000", "#0072B2", "#946853", "#D55E00"))+
      scale_fill_manual(values=c("#0072B2", "#D55E00")) 
plot(p)
ggsave("structuredSIR_results.png", p, width=6, height=4)


# plot the adjusted/weighted prevalence for constant vs. the constant reassortment rate
p = ggplot(const.est.data)+
  geom_point(aes(x=weightedPrevalence, y=prevalence))+
  geom_errorbar(aes(x=weightedPrevalence, ymin=prevalencel, ymax=prevalenceu), alpha=0.5) +
  geom_errorbar(aes(x=weightedPrevalence, ymin=prevalencel, ymax=prevalenceu), alpha=0.5) +
  # geom_text(aes(x=weightedPrevalence, y=prevalence, label=run))+
  theme_minimal() +
  # coord_cartesian(ylim=c(0,0.2)) + 
  geom_abline(intercept=0, slope=1, linetype="dashed")
  # scale the state colors only using alterating black and white
plot(p)



# plot the reassortment events inferred for each run and method
# sorted by the true number of reassortment events, also plot those
# true values, use one x axis slot for each run
p = ggplot(reassortment.data[order(reassortment.data$true),], aes(x=true, y=events))+
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      geom_point(alpha=0.9) +
      geom_errorbar(aes(ymin=eventsl, ymax=eventsu), alpha=0.5) +
      theme_minimal() +
      facet_wrap(method~., ncol=1)
plot(p)