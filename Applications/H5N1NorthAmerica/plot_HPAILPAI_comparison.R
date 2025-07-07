library(stringr)
library(seqinr)
library(ggplot2)
library("colorblindr")
library(ggtree)
library(treeio)
library(ggnewscale)
library(cowplot)

# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# plot the rate over time
reassortment_rate = data.frame()
# define the mrsi
mrsi = as.Date(c("2025-02-17", "2024-08-22", "2025-02-17", "2024-08-22"))

# define the segment order
segment_order = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

rate_shift_str = '0 0.12 0.24 0.36 0.48 0.6 0.72 0.84 0.96 1.08 1.2 1.32 1.44 1.56 1.68 1.8 1.92 2.04 2.16 2.28 2.4 2.52 2.64 2.76 2.88 3 3.12 3.24 3.36 3.48 3.6 3.72 3.84 3.96 4.08 4.2 4.32 4.44 4.56 4.68 4.8 4.92 5.04 5.16 5.28 5.4 5.52 5.64 5.76 5.88 6 6.8 7.6 8.4 9.2 10 1000 5000'
rate_shifts = as.numeric(strsplit(rate_shift_str, " ")[[1]])


lineage_colors <- c(
  HPAI    = "#E41A1C",  # red
  LPAI    = "#377EB8",  # blue
  unknown  = "#4DAF4A"  # green
)

clades = c("HPAI", "LPAI", "HPAI", "LPAI")
dataset = c("3 segments", "3 segments", "8 segments", "8 segments")
rerun = T


data = data.frame();
for (isIndependent in c(TRUE, FALSE)){
  

  if (isIndependent){
    if (rerun){
      # run log combined on the independent trees
      for (clade in clades){
        # run log combiner on the independent trees
        # system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
        #                "-burnin 50 -log ./out3seg/", clade, "_HLHxNx.independent.rep*.trees -o ./combined3seg/", clade, "_HLHxNx.independent.trees"))
        system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                       "-burnin 50 -log ./out3seg/", clade, "_HLHxNx.independent.rep*.log -o ./combined3seg/", clade, "_HLHxNx.independent.log"))
      }
    }
    # read in the log file for HLHxNx
    hpai_log_file <- read.csv("./combined3seg/HPAI_HLHxNx.independent.log", sep="\t")
    lpai_log_file <- read.csv("./combined3seg/LPAI_HLHxNx.independent.log", sep="\t")
    hpai_log_file_all_seg <- read.csv("./combined/HPAI_HLHxNx.independent.log", sep="\t")
    lpai_log_file_all_seg <- read.csv("./combined/LPAI_HLHxNx.independent.log", sep="\t")
    
  }else{
    
    if (rerun){
      for (clade in clades){
        # run log combiner on the dependent trees
        # system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
        #                "-burnin 50 -log ./out3seg/", clade, "_HLHxNx.dependent.rep*.trees -o ./combined3seg/", clade, "_HLHxNx.dependent.trees"))
        system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                       "-burnin 50 -log ./out3seg/", clade, "_HLHxNx.dependent.rep*.log -o ./combined3seg/", clade, "_HLHxNx.dependent.log"))
      }
    }
    # log_file <- read.csv("./combined3seg/HLHxNx.dependent.log", sep="\t")
    hpai_log_file <- read.csv("./combined3seg/HPAI_HLHxNx.dependent.log", sep="\t")
    lpai_log_file <- read.csv("./combined3seg/LPAI_HLHxNx.dependent.log", sep="\t")
    hpai_log_file_all_seg <- read.csv("./combined/HPAI_HLHxNx.dependent.log", sep="\t")
    lpai_log_file_all_seg <- read.csv("./combined/LPAI_HLHxNx.dependent.log", sep="\t")
  }
  c = 1
  for (log_file in list(hpai_log_file, lpai_log_file, hpai_log_file_all_seg, lpai_log_file_all_seg)){
    print(c)
    # get the clade from the log file name
    for (i in seq(1, length(rate_shifts))){
      # get the reassortment rate at this time point
      rate = log_file[, paste0("InfectedToRho.", i)]
      ne = log_file[, paste0("logNe.", i)]
      
      if (!isIndependent){
        rate = rate+ne
      }
      
      # get all the quantile from 0.05 to 1.0
      for (q in seq(0.05, 1.0, 0.05)){
        upper = quantile(rate, 1-q/2)
        lower = quantile(rate, q/2)
        if (q==1){
          lower = lower+0.03
          upper = upper-0.03
        }
        reassortment_rate = rbind(reassortment_rate, data.frame(
          time = mrsi[c]-rate_shifts[i]*365,
          quantile = q,
          upper = upper,
          lower = lower,
          isIndependent = isIndependent,
          clade = clades[c],
          dataset = dataset[c],
          rate = "reassortment"
        ))
        upper = quantile(ne, 1-q/2)
        lower = quantile(ne, q/2)
        if (q==1){
          lower = lower+0.03
          upper = upper-0.03
        }
        
        reassortment_rate = rbind(reassortment_rate, data.frame(
          time = mrsi[c]-rate_shifts[i]*365,
          quantile = q,
          upper = upper,
          lower = lower,
          isIndependent = isIndependent,
          clade = clades[c],
          dataset = dataset[c],
          rate = "Ne"
        ))
        
      }
    }
    
    clade = clades[c]
    c = c+1
    
  }
}


y_val = log(3)
y_off = 0.5

reassortment_rate$alpha = 0.2
# for q==1 set it to 1
reassortment_rate$alpha[reassortment_rate$quantile == 1] = 1.0



for (ind in c(TRUE, FALSE)){
  # plot the reassortment rate over time
  p = ggplot(reassortment_rate[reassortment_rate$isIndependent == ind & reassortment_rate$rate == "reassortment", ], 
             aes(x=time, y=upper, group=interaction(quantile, dataset), fill=dataset)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha=alpha)) +
    coord_cartesian(xlim=c(as.Date("2021-09-01"), max(mrsi)), ylim=c(-4, 2)) +
    scale_alpha(guide=F) +
    facet_wrap(.~clade, ncol = 1) +
    # scale_fill_manual(values=methods_colors, name="Method") +
    # scale_fill_manual(values=lineage_colors, name="Clade") +
    # mark the minimum and maximum time for each independent
    # label the y axis as exp(y)
    scale_y_continuous(breaks=c(log(0.05), log(0.1), log(0.2), log(0.4), log(0.8), log(1.6), log(3.2)), 
                       labels=c("0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "3.2")) +
    # scale_fill_OkabeIto()+
    xlab("Time") +
    ylab("Reassortment rate") +
    theme_minimal() 
  plot(p)
  if (ind){
    ggsave(p, filename="../../Figures/h5n1_hpai_lpai_rate_independent.pdf", width=6, height=2.5)
  }else{
    ggsave(p, filename="../../Figures/h5n1_hpai_lpai_reassortment_rate_dependent.pdf", width=6, height=2.5)
  }
  
  
  # read in positive cases
  cases = read.csv("./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv")
  cases$date = as.Date(cases$Date_Collected, format="%Y-%m-%d")
  
  # loop over all days
  min_date = min(cases$date)
  max_date = max(cases$date)
  smoothed_case_data = data.frame()
  for (d in seq(min_date+14, max_date-14, by="day")){
    # get all instances within that time window
    window = cases[cases$date >= d-7 & cases$date <= d+7, ]
    # get how many instances of Final_IAV are Detected
    total_AIV = sum(window$Final_IAV == "Detected", na.rm=TRUE)
    # get the number of high path cases 
    total_HPAI = sum(window$Final_Pathogenicity == "High Path AI", na.rm=TRUE)
    smoothed_case_data = rbind(smoothed_case_data, data.frame(
      date = d,
      positivity = (total_AIV-total_HPAI)/nrow(window),
      type = "LPAI"
    ))
    smoothed_case_data = rbind(smoothed_case_data, data.frame(
      date = d,
      positivity = total_HPAI/nrow(window),
      type = "HPAI"
    ))
  }
  smoothed_case_data$date = as.Date(smoothed_case_data$date, format="%Y-%m-%d")
  
  
  # plot the Ne over time, for both lineages
  p = ggplot(reassortment_rate[reassortment_rate$isIndependent == ind & reassortment_rate$rate == "Ne", ], 
             aes(x=time, y=upper, group=interaction(clade, quantile), fill=clade)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha=alpha)) +
    coord_cartesian(xlim=c(as.Date("2021-09-01"), max(mrsi))) +
    scale_alpha(guide=F) +
    # scale_fill_manual(values=methods_colors, name="Method") +
    scale_fill_manual(values=lineage_colors, name="Clade") +
    geom_line(data=smoothed_case_data, aes(x=date, y=positivity*12-1, color=type, group=type), method=NA, size=0.5) +
    scale_color_manual(values=c("HPAI"="#E41A1C", "LPAI"="#377EB8"), name="Type") +
    facet_wrap(.~dataset, ncol = 1) +
    xlab("Time") +
    ylab("Effective population size") +
    theme_minimal()
  plot(p)
  ggsave(p, filename=paste0("../../Figures/h5n1_hpai_lpai_ne_", ifelse(ind, "independent", "dependent"), ".pdf"), width=6, height=2.5)
  
}



