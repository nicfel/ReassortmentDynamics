library(stringr)
library(coda)
library(ggplot2)
# add library for color palette
library(RColorBrewer)
library("colorblindr")


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
CACHE_FILE = "SIR_plot_data.RDS"

if (file.exists(CACHE_FILE)) {
  cat("Loading cached data from", CACHE_FILE, "\n")
  cached = readRDS(CACHE_FILE)
  est.data = cached$est.data
  const.est.data = cached$const.est.data
  reassortment.data = cached$reassortment.data
  true.data = cached$true.data
  use.runs.1 = cached$use.runs.1
  use.runs.2 = cached$use.runs.2
} else {
# ---- DATA LOADING (skip when cache exists) ----
# date frame to save all results
est.data = data.frame()
const.est.data = data.frame()
# data frame to save the reassortment events
reassortment.data = data.frame()



# read in the file SIR_simulations.txt that contains the simulated values
simulated = read.table("SIR_simulations.txt", header=TRUE, sep="\t")

# get all .log files in /out
logfiles <- list.files("../out", pattern="*.log", full.names=TRUE)
# remove any files that contains structured
logfiles = logfiles[!grepl("structured", logfiles)]

# loop over all files
for (i in seq(1, length(logfiles))){
  t = try(read.table(logfiles[i], header=TRUE, sep="\t"))
  # if t is not a data.frame, skip the file
  if (!is.data.frame(t) || length(t$Sample) <10){
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
  
  # read in the corresponding xml by replacing .log with .xml and the folder with xmls
  xmlfile = gsub("out", "xmls", gsub(".log", ".xml", logfiles[i]))
  # read in the xml file
  xml = readLines(xmlfile)
  # look for the line with id="rateShifts"
  rateShifts = xml[grep("id=\"rateShifts\"", xml)]
  
  # split the line on " and get the second group
  rateShifts = as.numeric(strsplit(strsplit(rateShifts, "\"")[[1]][6], split=" ")[[1]])
  # time_points2 = seq(0,max(rateShifts), length.out=500)
  print(rateShifts)
  # split on . and get the second group as the method used
  method = strsplit(logfiles[i], "\\.")[[1]][4]
  print(method)
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][5])
  # if the method is not constant
  if (grepl("skygrowthNe", method)){
    transmision0 = simulated[run, "transmission"]
    
    for (q in seq(0,0.95, 0.05)){
      # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
      prevalencel = c()
      prevalenceu = c()
      time = c()
      
      # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
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
                                            quantile = q))     # get the value in simulated in the run row and "transmission" column
    }
  }  else if (grepl("skygrowth", method)){
    transmision0 = simulated[run, "transmission"]
    
    for (q in seq(0,0.95, 0.05)){
      # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
      prevalencel = c()
      prevalenceu = c()
      time = c()
      
      # loop over all the labels in t that start with reassortment%d, where %d is a number in time_points2
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
                                            quantile = q))     # get the value in simulated in the run row and "transmission" column
    }
    
  }else{
    # if the method is constant
    # get the value in simulated in the run row and "transmission" column
    transmision0 = simulated[run, "transmission"]
    # get the reassortment Rate
    prevalence = median(t[, "reassortmentRate"])
    # get the 2.5 and 97.5 quantile
    prevalencel = quantile(t[, "reassortmentRate"], 0.025)
    prevalenceu = quantile(t[, "reassortmentRate"], 0.975)
    # take the 25 and 75 quantile
    prevalencel5 = quantile(t[, "reassortmentRate"], 0.25)
    prevalenceu5 = quantile(t[, "reassortmentRate"], 0.75)
    prevalencel8 = quantile(t[, "reassortmentRate"], 0.9)
    prevalenceu8 = quantile(t[, "reassortmentRate"], 0.1)
    
    # dubblicate the above vectors to have the same value two times
    prevalence = c(prevalence, prevalence)
    prevalencel = c(prevalencel, prevalencel)
    prevalenceu = c(prevalenceu, prevalenceu)
    prevalencel5 = c(prevalencel5, prevalencel5)
    prevalenceu5 = c(prevalenceu5, prevalenceu5)
    # time
    time = c(0, 16)
    # add to the data.frame
    const.est.data = rbind(const.est.data, data.frame(method=method, run=run, time=time, 
            prevalence=prevalence/transmision0, prevalencel=prevalencel/transmision0, 
            prevalenceu=prevalenceu/transmision0,
            prevalencel5=prevalencel5/transmision0, prevalenceu5=prevalenceu5/transmision0,
            prevalencel8=prevalencel8/transmision0, prevalenceu8=prevalenceu8/transmision0))
  }
  reassortment.data = rbind(reassortment.data, data.frame(method=method, run=run, 
          events=median(t[, "network.reassortmentNodeCount"]), 
          eventsl=quantile(t[, "network.reassortmentNodeCount"], 0.025), 
          eventsu=quantile(t[, "network.reassortmentNodeCount"], 0.975)))
  
}
est.data$mrsi = NA
reassortment.data$true = NA
print("read in SIR files")


# get the true values for the prevalence by readin in the files in /master
true.data = data.frame()
# get all .log files in .master
logfiles <- list.files("./master", pattern="*.log", full.names=TRUE)

# loop over all files
for (i in seq(1, length(logfiles))){
# for (i in seq(1, 10)){
  #
  # split on . and _ simultanously and get the third group as the runnumber
  run = as.numeric(strsplit(logfiles[i], "\\.|_")[[1]][4])
  # if run is not part of unique(reassortment.data$run) skip
  if (!(run %in% unique(reassortment.data$run))){
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
  
  # compute the number of lineages over time by making a data frame
  currlins = 1 # tracks the number of lineages at any given time
  currtime = as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[2]]) # keeps the current time
  lineages_frame = data.frame() # dataframe for calcs after
  for (j in seq(length(lins),1,-1)){
    # split the string on : and get the second group
    tmp = strsplit(lins[[j]], split=":")[[1]]
    # add to the lineages data frame
    lineages_frame = rbind(lineages_frame, data.frame(start=currtime, end=as.numeric(tmp[[2]]), lineages=currlins))
    currtime = as.numeric(tmp[[2]])
    if (tmp[[1]] == "2"){
      currlins = currlins + 1
    }else{
      currlins = currlins - 1
    }
  }
  
  # get the most recent sampled individual
  mrsi = strsplit(lins[[1]], split=":")[[1]][[2]]
  # add the mrsi to the est.data
  est.data[est.data$run==run, "mrsi"] = as.numeric(mrsi)
  
  # Compute the number of infected individuals over time
  I <- 1
  time_I <- 0
  # Keep also track of the co-infection event times
  co_inf_times <- numeric()

  for (j in seq_along(sir)) {
      tmp <- str_split(sir[j], ':')[[1]]
      if (tmp[1] == '0') {
          I <- c(I, I[length(I)])
          time_I <- c(time_I, as.numeric(tmp[2]))
          co_inf_times <- c(co_inf_times, time_I[length(time_I)])
      } else if (tmp[1] == '1') {
          I <- c(I, I[length(I)] + 1)
          time_I <- c(time_I, as.numeric(tmp[2]))
      } else if (tmp[1] == '2') {
          I <- c(I, I[length(I)] - 1)
          time_I <- c(time_I, as.numeric(tmp[2]))
      } else if (tmp[1] == '3') {
          I <- c(I, I[length(I)])
          time_I <- c(time_I, as.numeric(tmp[2]))
      }
  }

  # for each start and end time in lineages_frame, compute the prevalence, waited by the number of lineages and the time.
  # keep track of the total time*lineages to normalize in the end
  normalization_factor = 0
  tot_prevalence = 0
  for (j in seq(1, length(lineages_frame$start))){
    # get all elements in time_I that are between start and end
    tmp = time_I[time_I >= lineages_frame$start[j] & time_I <= lineages_frame$end[j]]
    # get the average I's multiplied by the time intervals (diff in time_I)
    prevalence = sum(I[time_I >= lineages_frame$start[j] & time_I < lineages_frame$end[j]] * diff(tmp))
    # multiply by the number of lineages
    tot_prevalence = tot_prevalence + prevalence * lineages_frame$lineages[j]
    # compute the normalization factor
    normalization_factor = normalization_factor + lineages_frame$lineages[j] * sum(diff(tmp))
  }

  # get the true number of reassortment events in this run by counting
  # how often the strsplit(lins, split=":")[[:]][[1]] starts with 1:
  reassortment.data[reassortment.data$run==run, "true"] = sum(sapply(strsplit(lins, split=":"), function(x) x[[1]])=="1")

  # Store the results in lists
  I_list[[i]] <- I
  time_I_list[[i]] <- time_I
  co_inf_times_list[[i]] <- co_inf_times

  # get the population size
  popSize = simulated[run, "population_size"]
  
  # take ~100 points evenly distributed across time (not by index)
  target_times = seq(min(time_I), max(time_I), length.out=100)
  index = sort(unique(sapply(target_times, function(t) which.min(abs(time_I - t)))))
  
  # add the data to the data frame
  true.data = rbind(true.data, data.frame(run=run, time=time_I[index], prevalence=I[index], popSize=popSize, experienced_prev=tot_prevalence/normalization_factor))
}

# ggplot(lineages_frame, aes(x=start, xend=end, y=lineages, yend=lineages))+geom_segment()
dynamic_est_data = est.data[est.data$method!="constant", ]
# get all values of run in true.data$run
uni.run = unique(dynamic_est_data$run)
# for each value in uni.run, check if two methods are present, otherwise remove the value from uni.run
use.runs = c()
for (i in uni.run){
  if (length(unique(dynamic_est_data[dynamic_est_data$run==i,]$method))==2){
    use.runs = c(use.runs, i)
  }
}

# set.seed(15234)
# # pick 12 random runs
use.runs.1 = sample(use.runs[use.runs<=10], 5)
use.runs.2 = sample(use.runs[use.runs>10], 5)

  # Save for quick plotting iteration
  saveRDS(list(est.data=est.data, const.est.data=const.est.data, reassortment.data=reassortment.data,
               true.data=true.data, use.runs.1=use.runs.1, use.runs.2=use.runs.2), CACHE_FILE)
  cat("Saved data to", CACHE_FILE, "\n")
}

# ---- PLOTTING ----
dynamic_est_data = est.data[est.data$method!="constant", ]

# plot the results
dd1 = dynamic_est_data[dynamic_est_data$run %in% use.runs.1, ]
xrange1 = range(dd1$mrsi - dd1$time, na.rm=TRUE)
xpad1 = max(diff(xrange1) * 0.02, 0.1)
p = ggplot(true.data[true.data$run %in% use.runs.1,], aes(x=time, y=prevalence/popSize))+
      geom_line() +
      geom_ribbon(data=dd1, aes(x=mrsi-time, ymin=pmin(prevalencel, 1), ymax=pmin(prevalenceu, 1), fill=method, group=quantile), alpha=0.1, inherit.aes=FALSE) +
      facet_wrap(method~run, ncol=5, scales='free') +
      theme_minimal() +
      theme(strip.text=element_blank()) +
      coord_cartesian(xlim=c(xrange1[1]-xpad1, xrange1[2]+xpad1)) +
      scale_color_manual(values=c("#0072B2", "#D55E00"))+
      scale_fill_manual(values=c("#0072B2", "#D55E00")) +
      guides(fill=guide_legend(override.aes=list(alpha=1))) +
      xlab("Time") +
      ylab("Prevalence")
plot(p)
# save the results
ggsave("./../../Figures/Sir.pdf", p, width = 9, height = 6)

dd2 = dynamic_est_data[dynamic_est_data$run %in% use.runs.2, ]
xrange2 = range(dd2$mrsi - dd2$time, na.rm=TRUE)
xpad2 = max(diff(xrange2) * 0.02, 0.1)
p = ggplot(true.data[true.data$run %in% use.runs.2 ,], aes(x=time, y=prevalence/popSize))+
  geom_line() +
  geom_ribbon(data=dd2, aes(x=mrsi-time, ymin=pmin(prevalencel, 1), ymax=pmin(prevalenceu, 1), fill=method, group=quantile), alpha=0.1, inherit.aes=FALSE) +
  facet_wrap(method~run, ncol=5, scales='free') +
  theme_minimal() +
  theme(strip.text=element_blank()) +
  coord_cartesian(xlim=c(xrange2[1]-xpad2, xrange2[2]+xpad2)) +
  scale_color_manual(values=c("#0072B2", "#D55E00"))+
  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
  guides(fill=guide_legend(override.aes=list(alpha=1))) +
  xlab("Time") +
  ylab("Prevalence")
plot(p)
# save the results
ggsave("./../../Figures/Sir_superspreading.pdf", p, width = 9, height = 6)


# plot the constant results by first making a dataframe with the constant results form true.data and est.data by matching the run
const_data = data.frame()
for (i in unique(true.data$run)){
  # get the ratio between true and estimated number of events
  if (length(reassortment.data[reassortment.data$run==i & reassortment.data$method=="constant",]$events)==0){
    next
  }
  ratio = reassortment.data[reassortment.data$run==i & reassortment.data$method=="constant",]
  
  # get the corresponding row in true.data that has method constant 
  # and the run i, but only pick the first row for which this is true
  const_data = rbind(const_data, data.frame(run=i, experienced_prev=true.data[true.data$run==i,]$experienced_prev[1],
          popSize=true.data[true.data$run==i,]$popSize[1], 
          prevalence=const.est.data[const.est.data$run==i,]$prevalence[1],
          prevalencel=const.est.data[const.est.data$run==i,]$prevalencel[1],
          prevalenceu=const.est.data[const.est.data$run==i,]$prevalenceu[1],
          ratio = ratio$events/ratio$true))  
}

const_data$super = "with"
const_data[const_data$run<=10, "super"] = "without"
# make an x equal y plot
p = ggplot(const_data[order(const_data$experienced_prev/const_data$popSize),], 
           aes(x=experienced_prev/popSize*ratio, y=prevalence, color=super))+
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      geom_point(alpha=0.9) +
      geom_errorbar(aes(ymin=prevalencel, ymax=prevalenceu), alpha=0.5) +
  xlab("True, time weighted prevalence")+
  ylab("Estimated prevalence from reassortment rate")+
  scale_color_OkabeIto(name="superspreading")+
      theme_minimal()
plot(p)
ggsave("./../../Figures/Sir_constant.pdf", p, width = 6, height = 3)


# plot the reassortment events inferred for each run and method
# sorted by the true number of reassortment events, also plot those
# true values, use one x axis slot for each run
p = ggplot(reassortment.data[order(reassortment.data$true),], aes(x=true, y=events))+
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      geom_point(alpha=0.9) +
      geom_errorbar(aes(ymin=eventsl, ymax=eventsu), alpha=0.5) +
      theme_minimal() +
      xlab("True number of reassortment events")+
      ylab("Estimated number of reassortment events")+
      facet_wrap(method~., ncol=1)
plot(p)
ggsave("./../../Figures/Sir_constant_events.pdf", p, width = 6, height = 6)


p = ggplot(reassortment.data, aes(x=events-true))+
  geom_vline(xintercept=0,linetype="dashed") +
  geom_histogram(alpha=0.9, binwidth = 5) +
  theme_minimal() +
  facet_wrap(method~., ncol=1) +
  xlab("estimated - true number of reassortment events")
plot(p)
ggsave("./../../Figures/Sir_constant_events_distr.pdf", p, width = 6, height = 6)

