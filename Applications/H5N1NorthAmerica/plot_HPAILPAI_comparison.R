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
mrsi = as.Date(c("2025-02-17", "2024-08-14"))

# define the segment order
segment_order = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

rate_shift_str = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5 10 15 20 25 30 300'
rate_shifts = as.numeric(strsplit(rate_shift_str, " ")[[1]])


lineage_colors <- c(
  HPAI    = "#E41A1C",  # red
  LPAI    = "#377EB8",  # blue
  unknown  = "#4DAF4A"  # green
)

clades = c("HPAI", "LPAI")
rerun = F


data = data.frame();
for (isIndependent in c(TRUE, FALSE)){
  

  if (isIndependent){
    if (rerun){
      # run log combined on the independent trees
      for (clade in clades){
        # run log combiner on the independent trees
        # system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
        #                "-burnin 50 -log ./out3seg/", clade, "_450.independent.rep*.trees -o ./combined3seg/", clade, "_450.independent.trees"))
        system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                       "-burnin 50 -log ./out3seg/", clade, "_450.independent.rep*.log -o ./combined3seg/", clade, "_450.independent.log"))
      }
    }
    # read in the log file for 450
    hpai_log_file <- read.csv("./combined3seg/HPAI_450.independent.log", sep="\t")
    lpai_log_file <- read.csv("./combined3seg/LPAI_450.independent.log", sep="\t")
  }else{
    
    if (rerun){
      for (clade in clades){
        # run log combiner on the dependent trees
        # system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
        #                "-burnin 50 -log ./out3seg/", clade, "_450.dependent.rep*.trees -o ./combined3seg/", clade, "_450.dependent.trees"))
        system(paste0("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                       "-burnin 50 -log ./out3seg/", clade, "_450.dependent.rep*.log -o ./combined3seg/", clade, "_450.dependent.log"))
      }
    }
    # log_file <- read.csv("./combined3seg/450.dependent.log", sep="\t")
    hpai_log_file <- read.csv("./combined3seg/HPAI_450.dependent.log", sep="\t")
    lpai_log_file <- read.csv("./combined3seg/LPAI_450.dependent.log", sep="\t")
  }
  c = 1
  for (log_file in list(hpai_log_file, lpai_log_file)){
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
  p = ggplot(reassortment_rate[reassortment_rate$isIndependent == ind & reassortment_rate$rate == "reassortment", ], aes(x=time, y=upper, group=quantile, fill=clade)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha=alpha), fill="grey29") +
    coord_cartesian(xlim=c(as.Date("2021-09-01"), max(mrsi)), ylim=c(-4, 2)) +
    scale_alpha(guide=F) +
    facet_wrap(.~clade, ncol = 1) +
    # scale_fill_manual(values=methods_colors, name="Method") +
    scale_fill_manual(values=lineage_colors, name="Clade") +
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
    
    xlab("Time") +
    ylab("Effective population size") +
    theme_minimal()
  plot(p)
  ggsave(p, filename=paste0("../../Figures/h5n1_hpai_lpai_ne_", ifelse(ind, "independent", "dependent"), ".pdf"), width=6, height=2.5)
  
}



