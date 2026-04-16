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

# Set the directory to the directory of the file (works with Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
match <- grep(file_arg, args)
if (length(match) > 0) {
  this.dir <- dirname(normalizePath(sub(file_arg, "", args[match])))
} else {
  this.dir <- getwd()
}
setwd(this.dir)

# Cache file for quick plotting iteration (delete to regenerate from source data)
CACHE_FILE = "StructuredSIR_plot_data.RDS"

if (file.exists(CACHE_FILE)) {
  cat("Loading cached data from", CACHE_FILE, "\n")
  cached = readRDS(CACHE_FILE)
  est.data = cached$est.data
  const.est.data = cached$const.est.data
  reassortment.data = cached$reassortment.data
  true.data = cached$true.data
  use.runs.1 = cached$use.runs.1
  reassortment.events = cached$reassortment.events
} else {

# date frame to save all results
est.data = data.frame()
const.est.data = data.frame()
# data frame to save the reassortment events
reassortment.data = data.frame()

#define the number of states
nstates = 50

# read in the file SIR_simulations.txt that contains the simulated values
simulated = read.table("structuredSIR_simulations.txt", header=TRUE, sep="\t")

# get all .log files in /out
logfiles <- list.files("../out", pattern="*.log", full.names=TRUE)
logfiles <- logfiles[grepl("structured", logfiles)]


# remove all files that do not contain the word variable not variable
# logfiles = logfiles[grepl("variable", logfiles) ]

# remove all log files that are not run 13
# logfiles = logfiles[grepl("_13.", logfiles) ]

# loop over all files
for (i in seq(1, length(logfiles))){
  # read in the log file to get the network height
  filename = gsub("../out/", "./master/", logfiles[[i]])  
  filename = gsub(".skygrowthNe.", ".", filename)  
  filename = gsub(".skygrowth.", ".",filename)  
  filename = gsub(".constant.", ".", filename)  

  lines = readLines(filename)
  data = strsplit(lines[[2]], split="\t")[[1]]
  lins <- str_split(str_replace_all(data[[3]], "\\[|\\]", ""), ', ')[[1]]
  
  timediff = as.numeric(strsplit(lins[[1]], split=":")[[1]][[2]]) - 
    as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[2]])

  t = read.table(logfiles[i], header=TRUE, sep="\t")
  # remove 10 % as burnin
  t = t[round(0.2*nrow(t)):nrow(t),]
  # calculate the ESS for columns 2 to 5
  ess = sapply(t[,2:5], function(x) effectiveSize(x))
  # calculate the minimum ESS
  minESS = min(ess)
  # if the minimum ESS is below 100, print the filename
  # if the minimum ESS is below 100, print the filename
  if (minESS < 5) {
    print(minESS)
    print(logfiles[i])
    next
  }
  
  # read in the corresponding xml by replacing .log with .xml and the folder with xmls
  xmlfile = gsub("out", "xmls", gsub(".log", ".xml", logfiles[i]))
  # read in the xml file
  xml = readLines(xmlfile)
  # look for the line with id="rateShifts"
  rateShifts = xml[grep("id=\"rateShifts\"", xml)]
  
  # split the line on " and get the second group
  rateShifts = as.numeric(strsplit(strsplit(rateShifts, "\"")[[1]][6], split=" ")[[1]])
  time_points2 = seq(0,max(rateShifts), length.out=500)
  
  # split on . and get the second group as the method used
  method = strsplit(logfiles[i], "\\.")[[1]][4]
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][5])
  # if the method is not constant
  if (grepl("skygrowthNe", method)){
    # get the value in simulated in the run row and "transmission" column
    transmision0 = simulated[run, "transmission"]
    # initialize the prevalence vector
    time = c()
    
    for (q in seq(0,0.95, 0.05)){
      # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
      prevalencel = c()
      prevalenceu = c()
      for (j in seq(1, length(rateShifts), 1)){
        vals = t[,paste0("logNe.", j)] + t[,paste0("InfectedToRho.", j)]
        prevalencel = c(prevalencel, quantile(vals, (1-q)/2))
        prevalenceu = c(prevalenceu, quantile(vals, (1+q)/2))
        # time
        time = c(time, rateShifts[j])
      }
      # add the data to a data frame
      est.data = rbind(est.data, data.frame(method=method, run=run, time=time, 
                                            prevalencel= exp(prevalencel)/transmision0, prevalenceu=exp(prevalenceu)/transmision0,
                                            quantile = q)) 
    }
    
  }  else if (grepl("skygrowth", method)){
      # get the value in simulated in the run row and "transmission" column
      transmision0 = simulated[run, "transmission"]
      # initialize the prevalence vector
      time = c()
      
      for (q in seq(0,0.95, 0.05)){
        # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
        prevalencel = c()
        prevalenceu = c()
        for (j in seq(1, length(rateShifts), 1)){
          vals = t[,paste0("InfectedToRho.", j)]
          prevalencel = c(prevalencel, quantile(vals, (1-q)/2))
          prevalenceu = c(prevalenceu, quantile(vals, (1+q)/2))
          # time
          time = c(time, rateShifts[j])
        }
        # add the data to a data frame
        est.data = rbind(est.data, data.frame(method=method, run=run, time=time, 
                                              prevalencel= exp(prevalencel)/transmision0, prevalenceu=exp(prevalenceu)/transmision0,
                                              quantile = q)) 
      }
    }else{
    # if the method is constant
    # get the value in simulated in the run row and "transmission" column
    transmision0 = simulated[run, "transmission"]
    time = c(0, 16)
    
    # get the reassortment Rate
    for (q in seq(0,0.95, 0.05)){
      vals = t[,"reassortmentRate"]
      # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
      prevalencel = quantile(vals, (1-q)/2)
      prevalenceu = quantile(vals, (1+q)/2)
      prevalencel = c(prevalencel, quantile(vals, (1-q)/2))
      prevalenceu = c(prevalenceu, quantile(vals, (1+q)/2))

      # add the data to a data frame
      const.est.data = rbind(const.est.data, data.frame(method=method, run=run, time=time, 
                                            prevalencel= prevalencel/transmision0, prevalenceu=prevalenceu/transmision0,
                                            quantile = q)) 
    }
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

# loop over all files
for (i in seq(1, length(logfiles))){
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
  
  print('lins done')
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

  time_I <- c()
  # Keep also track of the co-infection event times
  co_inf_times <- numeric()  
  
  # read in the xml and find the line with <populationSize spec="IntegerParameter" value="..."", save the populaiton sizes into a vector
  xmlfile = gsub(".log", ".xml", logfiles[i])
  xml = readLines(xmlfile)
  popSizeCurrent = as.numeric(strsplit(strsplit(xml[grep("<populationSize", xml)], split="\"")[[1]][4], split=" ")[[1]])
  popSize = c()

  vals <- str_split(sir, ':')
  # do this for the 4 group over everything in vals all at once
  type_vals = sapply(vals, function(x) as.numeric(strsplit(x[[4]], split="_")[[1]][[1]])+1)
  k=2
  first =TRUE
  tot_prevalence = 0
  normalization_factor = 0
  for (j in seq_along(sir)) {
    tmp <- vals[[j]]
    type = type_vals[j]
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
    t = as.numeric(tmp[2]) 
    if (t>= lineages_frame$start[k]){
      if (first){
        I = Icurrent
        popSize = popSizeCurrent
        time_I = c(t)
        first = FALSE
      }else{
        I = rbind(I, Icurrent)
        popSize = rbind(popSize, popSizeCurrent)
        time_I = c(time_I, t)
      }
    }
    if (t > lineages_frame$end[k]){
      #calculate teh ratio
      ratio = I / popSize
      lins_state <- as.matrix(lineages_frame[k, 3:ncol(lineages_frame)])
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
      
      total_time = time_I[length(time_I)] - time_I[1]
      
      
      if (total_time>0){
        time_weighted_prevalence = sum(weighted_prevalence * diff(time_I))/total_time
        tot_prevalence = tot_prevalence + time_weighted_prevalence * total_time * sum(lins_state)
        # compute the normalization factor to later normalize by the total time
        normalization_factor = normalization_factor + sum(lins_state) * sum(diff(time_I))
      }else{
        time_weighted_prevalence = mean(weighted_prevalence);
      }
      
      
      true.data = rbind(true.data, data.frame(run=run, time=t, 
                                              weighted_prev=time_weighted_prevalence,
                                              mean_prev=mean(ratio),
                                              max_prev= max(ratio),
                                              num_lineages=sum(lins_state),
                                              method="time varying"))
      
      
      # if the time is greater than the end time, go to the next line
      k = k + 1
      
      if (k>nrow(lineages_frame)){
        break
      }
      # reset I, popSize and time_I
      I = matrix(0, nrow=0, ncol=nstates)
      popSize = matrix(0, nrow=0, ncol=nstates)
      time_I = c()
    }
  }
  

  print('sir done')
  
  # # for each start and end time in lineages_frame, compute the prevalence, 
  # # waited by the number of lineages in a given location and the time.
  # # keep track of the total time*lineages to normalize in the end
  # normalization_factor = 0
  # tot_prevalence = 0
  # for (j in seq(2, length(lineages_frame$start))){
  #   # get all elements in time_I that are between start and end to get the true number of infected for that time period
  #   tmp = time_I[time_I >= lineages_frame$start[j] & time_I <= lineages_frame$end[j]]
  #   # multiply the number of infected individuals in each state by the number of lineages
  #   # in the same state for lineages_frame[j,3:end], but in R language
  #   infected_state = I[time_I >= lineages_frame$start[j] & time_I < lineages_frame$end[j],]
  #   # multiply the matrix infected_state with the dataframe vector lins_state,
  #   # but fix the orientation of lins_state
  #   popSizeState = popSize[time_I >= lineages_frame$start[j] & time_I < lineages_frame$end[j],]
  #   ratio = infected_state / popSizeState
  #   
  #   lins_state <- as.matrix(lineages_frame[j, 3:ncol(lineages_frame)])
  #   # divide the number of infected individuals in each state by the population size
  #   # in the same state
  #   weighted_prevalence = c()
  #   if (length(ratio)>nstates){
  #     for (row in seq(1, length(ratio)/nstates)){
  #       weighted_prevalence[row] = sum(ratio[row,] * lins_state)
  #     }
  #     weighted_prevalence = weighted_prevalence/sum(lins_state)
  #   }else if (length(ratio)==nstates){
  #     weighted_prevalence = sum(ratio * lins_state)
  #     weighted_prevalence = weighted_prevalence/sum(lins_state)
  #   }
  # 
  # 
  #   # check if any elemente in weighted_prevalence is below 0
  #   if (any(weighted_prevalence<0)){
  #     print("negative prevalence")
  #     dsa
  #   }
  # 
  #   # now also weigh for the time of each interval
  #   total_time = tmp[length(tmp)] - tmp[1]
  #   # weight the prevalence_times_lineages by the time of each interval
  #   if (total_time>0){
  #     time_weighted_prevalence = sum(weighted_prevalence * diff(tmp))/total_time
  #     tot_prevalence = tot_prevalence + time_weighted_prevalence * total_time * sum(lins_state)
  #     # compute the normalization factor to later normalize by the total time
  #     normalization_factor = normalization_factor + sum(lins_state) * sum(diff(tmp))
  #   }else{
  #     time_weighted_prevalence = mean(weighted_prevalence);
  #   }
  #   
  #   if (length(tmp)>1){
  #     true.data = rbind(true.data, data.frame(run=run, time=(tmp[length(tmp)]+tmp[1])/2, 
  #                                             weighted_prev=time_weighted_prevalence,
  #                                             mean_prev=mean(ratio),
  #                                             max_prev= max(ratio),
  #                                             num_lineages=sum(lins_state),
  #                                             method="time varying"))
  #   }
  #   else{
  #     true.data = rbind(true.data, data.frame(run=run, time=tmp[1], 
  #                                             weighted_prev=time_weighted_prevalence,
  #                                             mean_prev=mean(ratio),
  #                                             max_prev= max(ratio),
  #                                             num_lineages=sum(lins_state),
  #                                             method="time varying"))
  #   }
  # }
  # 
  print(tot_prevalence/normalization_factor)
  # add the weighted prevalence to the const.est.data$weightedPrevalence by matching runs
  const.est.data[const.est.data$run==run, "weightedPrevalence"] = tot_prevalence/normalization_factor

  # experienced prevalence for constant calculation experienced_prev=tot_prevalence/normalization_factor)
  # get the true number of reassortment events in this run by counting
  # how often the strsplit(lins, split=":")[[:]][[1]] starts with 1:
  reassortment.data[reassortment.data$run==run, "true"] = sum(sapply(strsplit(lins, split=":"), function(x) x[[1]])=="1")

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
  maxu = max(est.data[est.data$run==run & (est.data$method=="skyline"), "prevalenceu"])
  reassortment.events[reassortment.events$run==run, "y"] = maxu
}


uni.run = unique(est.data$run)
# for each value in uni.run, check if two methods are present, otherwise remove the value from uni.run
use.runs = uni.run
# for (i in uni.run){
#   if (length(unique(est.data[est.data$run==i,]$method))==2){
#     use.runs = c(use.runs, i)
#   }
# }
set.seed(15234)
# pick 12 random runs
use.runs.1 = use.runs

  # Save for quick plotting iteration
  saveRDS(list(est.data=est.data, const.est.data=const.est.data, reassortment.data=reassortment.data,
               true.data=true.data, use.runs.1=use.runs.1, reassortment.events=reassortment.events), CACHE_FILE)
  cat("Saved data to", CACHE_FILE, "\n")
}

# ---- PLOTTING ----
subet_data = data.frame()
# for each run in true data, select ~100 points evenly distributed across time
for (i in unique(true.data$run)){
  tp = true.data[true.data$run==i,]
  target_times = seq(min(tp$time), max(tp$time), length.out=100)
  index = sort(unique(sapply(target_times, function(t) which.min(abs(tp$time - t)))))
  subet_data = rbind(subet_data, tp[index,])
}

# plot the results, for true.data, plot the stacked prevalence over time
subet_data$method="skygrowth"
dd_sg = est.data[is.element(est.data$run, use.runs.1) & est.data$method=="skygrowth", ]
xrange_sg = range(dd_sg$mrsi - dd_sg$time, na.rm=TRUE)
xpad_sg = max(diff(xrange_sg) * 0.02, 0.1)
# reassortment events at max y per facet (run)
reassortment_sg = reassortment.events[is.element(reassortment.events$run, use.runs.1), c("run", "time")]
reassortment_sg$method = "skygrowth"
for (run in use.runs.1) {
  idx = reassortment_sg$run == run
  subet_max = max(subet_data[subet_data$run==run, c("weighted_prev", "mean_prev")], na.rm=TRUE)
  est_max = if (any(dd_sg$run==run)) max(pmin(dd_sg[dd_sg$run==run, "prevalencel"], 1), pmin(dd_sg[dd_sg$run==run, "prevalenceu"], 1), na.rm=TRUE) else 0
  reassortment_sg$y[idx] = max(subet_max, est_max, na.rm=TRUE)
}
p = ggplot(subet_data[is.element(subet_data$run, use.runs.1),])+
      geom_ribbon(data=dd_sg, aes(x=mrsi-time, ymin=pmin(prevalencel, 1), ymax=pmin(prevalenceu, 1), fill=method, group=quantile), alpha=0.10, inherit.aes=FALSE) +
      geom_line(aes(x=time, y=weighted_prev, color="lineage weighted mean prevalence across states", method="skygrowth"))+
      geom_line(aes(x=time, y=mean_prev, color="mean prevalence across states", method="skygrowth"))+
      facet_wrap(method~run, ncol=5, scales = "free") +
      theme_minimal() +
      theme(strip.text=element_blank(), legend.position="top") +
      coord_cartesian(xlim=c(xrange_sg[1]-xpad_sg, xrange_sg[2]+xpad_sg)) +
      guides(fill=guide_legend(override.aes=list(alpha=1))) +
      xlab("Time") +
      ylab("Prevalence")+
      geom_count(data=reassortment_sg, aes(x=time, y=y), alpha=0.1, size=2)+
      
      scale_color_manual(values=c("#946853", "#000000", "#0072B2", "#946853", "#D55E00"))+
      scale_fill_manual(values=c("#0072B2")) 
plot(p)
ggsave("./../../Figures/StructuredSir_skygrowth.pdf", p, width = 9, height = 5)


# plot the results, for true.data, plot the stacked prevalence over time
subet_data$method="skygrowthNe"
dd_sgne = est.data[is.element(est.data$run, use.runs.1) & est.data$method=="skygrowthNe", ]
xrange_sgne = range(dd_sgne$mrsi - dd_sgne$time, na.rm=TRUE)
xpad_sgne = max(diff(xrange_sgne) * 0.02, 0.1)
reassortment_sgne = reassortment.events[is.element(reassortment.events$run, use.runs.1), c("run", "time")]
reassortment_sgne$method = "skygrowthNe"
for (run in use.runs.1) {
  idx = reassortment_sgne$run == run
  subet_max = max(subet_data[subet_data$run==run, c("weighted_prev", "mean_prev")], na.rm=TRUE)
  est_max = if (any(dd_sgne$run==run)) max(pmin(dd_sgne[dd_sgne$run==run, "prevalencel"], 1), pmin(dd_sgne[dd_sgne$run==run, "prevalenceu"], 1), na.rm=TRUE) else 0
  reassortment_sgne$y[idx] = max(subet_max, est_max, na.rm=TRUE)
}
p = ggplot(subet_data[is.element(subet_data$run, use.runs.1),])+
  geom_ribbon(data=dd_sgne, aes(x=mrsi-time, ymin=pmin(prevalencel, 1), ymax=pmin(prevalenceu, 1), fill=method, group=quantile), alpha=0.10, inherit.aes=FALSE) +
  geom_line(aes(x=time, y=weighted_prev, color="lineage weighted mean prevalence across states", method="skygrowthNe"))+
  geom_line(aes(x=time, y=mean_prev, color="mean prevalence across states", method="skygrowthNe"))+
  facet_wrap(method~run, ncol=5, scales = "free") +
  theme_minimal() +
  theme(strip.text=element_blank(), legend.position="top") +
  coord_cartesian(xlim=c(xrange_sgne[1]-xpad_sgne, xrange_sgne[2]+xpad_sgne)) +
  guides(fill=guide_legend(override.aes=list(alpha=1))) +
  xlab("Time") +
  ylab("Prevalence")+
  geom_count(data=reassortment_sgne, aes(x=time, y=y), alpha=0.2, size=2)+
  
  scale_color_manual(values=c("#946853", "#000000", "#0072B2", "#946853", "#D55E00"))+
  scale_fill_manual(values=c("#D55E00")) 
plot(p)
ggsave("./../../Figures/StructuredSir_skygrowthNe.pdf", p, width = 9, height = 5)





# plot the adjusted/weighted prevalence for constant vs. the constant reassortment rate
p = ggplot(const.est.data[const.est.data$quantile==0.95,])+
  geom_errorbar(aes(x=weightedPrevalence, ymin=prevalencel, ymax=prevalenceu ), width=0) +
  geom_point(data=const.est.data[const.est.data$quantile==0, ], aes(x=weightedPrevalence, y=prevalencel)) +
    # geom_text(aes(x=weightedPrevalence, y=prevalencel, label=run))+
  theme_minimal() +
  xlab("True, time weighted prevalence")+
  ylab("Estimated prevalence from reassortment rate")+
  geom_abline(intercept=0, slope=1, linetype="dashed")
  # scale the state colors only using alterating black and white
plot(p)
ggsave("./../../Figures/StructuredSir_constant_prevalence.pdf", p, width = 6, height = 6)



# plot the reassortment events inferred for each run and method
# sorted by the true number of reassortment events, also plot those
# true values, use one x axis slot for each run
p = ggplot(reassortment.data[order(reassortment.data$true),], aes(x=true, y=events))+
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      geom_point(alpha=0.9) +
      geom_errorbar(aes(ymin=eventsl, ymax=eventsu), alpha=0.5, width=0) +
      theme_minimal() +
      xlab("True number of reassortment events")+
      ylab("Estimated number of reassortment events")+
      facet_wrap(method~., ncol=1)
plot(p)

ggsave("./../../Figures/StructuredSir_constant_events.pdf", p, width = 6, height = 6)
