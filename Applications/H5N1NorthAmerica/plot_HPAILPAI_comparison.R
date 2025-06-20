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

rate_shift_str = '0 0.270743639921722 0.541487279843444 0.812230919765166 1.08297455968689 1.35371819960861 1.62446183953033 1.89520547945205 2.16594911937378 2.4366927592955 2.70743639921722 2.97818003913894 3.24892367906067 3.51966731898239 3.79041095890411 9.03232876712329 14.2742465753425 19.5161643835616 24.7580821917808 30'
rate_shifts = as.numeric(strsplit(rate_shift_str, " ")[[1]])


lineage_colors <- c(
  HPAI    = "#E41A1C",  # red
  LPAI    = "#377EB8",  # blue
  unknown  = "#4DAF4A"  # green
)

clades = c("HPAI", "LPAI")
rerun = TRUE

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
    # get the clade from the log file name
    for (i in seq(1, length(rate_shifts))){
      # get the reassortment rate at this time point
      rate = log_file[, paste0("InfectedToRho.", i)]
      if (!isIndependent){
        rate = rate+log_file[, paste0("logNe.", i)]
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
      }
      for (i in seq(1, length(rate_shifts))){
        # get the reassortment rate at this time point
        rate = log_file[, paste0("logNe.", i)]
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
            rate = "Ne"
          ))
        }
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
  
  # plot the Ne over time, for both lineages
  p = ggplot(reassortment_rate[reassortment_rate$isIndependent == ind & reassortment_rate$rate == "Ne", ], 
             aes(x=time, y=upper, group=interaction(clade, quantile), fill=clade)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha=alpha)) +
    coord_cartesian(xlim=c(as.Date("2021-09-01"), max(mrsi))) +
    scale_alpha(guide=F) +
    # scale_fill_manual(values=methods_colors, name="Method") +
    scale_fill_manual(values=lineage_colors, name="Clade") +
    # mark the minimum and maximum time for each independent
    # label the y axis as exp(y)
    # scale_y_continuous(breaks=c(log(0.05), log(0.1), log(0.2), log(0.4), log(0.8), log(1.6), log(3.2)), 
    #                    labels=c("0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "3.2")) +
    # scale_fill_OkabeIto()+
    xlab("Time") +
    ylab("Effective population size") +
    theme_minimal()
  plot(p)
  ggsave(p, filename=paste0("../../Figures/h5n1_hpai_lpai_ne_", ifelse(ind, "independent", "dependent"), ".pdf"), width=6, height=2.5)
  
}



