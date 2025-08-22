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
mrsi = as.Date("2025-02-17")
mrsi_hpai = mrsi
mrsi_lpai = as.Date("2024-08-22")

# define the segment order
segment_order = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

rate_shift_str = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5 10 15 20 25 30 1000'
rate_shifts = as.numeric(strsplit(rate_shift_str, " ")[[1]])

# 1) Methods (3 distinct hues from Set1)
methods_colors <- c(
  constant    = "#E41A1C",  # red
  skyline     = "#377EB8",  # blue
  skyline_Ne  = "#4DAF4A"   # green
)


# 3) Species (2 contrasting hues from Dark2)
species_colors <- c(
  cow  = "#238B45",  # teal
  bird = "#D95F02"   # burnt orange
)

# 3) Species (2 contrasting hues from Dark2)
clade_colors <- c(
  "B3.13"  = "#238B45",  # teal
  "D1.1" = "#D95F02"   # burnt orange
)



lineage_colors <- c(
  HPAI    = "#E41A1C",  # red
  LPAI    = "#377EB8",  # blue
  unknown  = "#4DAF4A"  # green
)

clades = c("B3.13", "D1.1")
rerun = T


data_type = c("hpai", "lpai", "total", "overlap", "independent", "dependent")

data = data.frame();
for (data_type in data_type){
  if (rerun){
    # run log combined on the independent trees
    system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                 "-burnin 20 -log ./out3/HPAI_HLHxNx.", data_type,"*.log -o ./combined/HPAI_HLHxNx.", data_type, ".log", sep=""))

  }
  # read in the log file for HLHxNx
  log_file <- read.table(paste("./combined/HPAI_HLHxNx.", data_type, ".log", sep=""), header=TRUE, comment.char="#")
  
  # read in the corresponding xml file in the xmls folder and look for logNe
  xml_file <- readLines(paste("./xmls/HPAI_HLHxNx.", data_type, ".rep0.xml", sep=""))
  # lookg for a line with logStandardizedCases and get the value
  logNe_line <- grep("logStandardizedCases", xml_file, value=TRUE)
  if (length(logNe_line) == 0) {
    logNe_value = rep(0, length(rate_shifts))
  }else{
    # extract the value from the line
    logNe_value <- as.numeric(strsplit(strsplit(logNe_line, '"')[[1]][[6]], " ")[[1]])
  }

  
  # choose the appropriate mrsi vector
  mrsi_tmp <- mrsi
  for (i in seq(1, length(rate_shifts))){
    # get the reassortment rate at this time point
    rate = log_file[, paste0("InfectedToRho.", i)]
    # rate = rate+ logNe_value[i]
    # if (data_type == "dependent"){
    #   rate = rate + log_file[, paste0("logNe.", i)]
    # }

    # get all the quantile from 0.05 to 1.0
    for (q in seq(0.05, 1.0, 0.05)){
      upper = quantile(rate, 1-q/2)
      lower = quantile(rate, q/2)
      if (q==1){
        lower = lower+0.03
        upper = upper-0.03
      }
      reassortment_rate = rbind(reassortment_rate, data.frame(
        time = mrsi_tmp - rate_shifts[i]*365,
        quantile = q,
        upper = upper,
        lower = lower,
        name = data_type
      ))
    }
  }
}
reassortment_rate$time = as.Date(reassortment_rate$time)

reassortment_rate$alpha = 0.2
# for q==1 set it to 1
reassortment_rate$alpha[reassortment_rate$quantile == 1] = 1.0

ggplot(reassortment_rate, aes(x=time, y=upper, group=interaction(name, quantile), fill=name, alpha=quantile)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  coord_cartesian(xlim=c(as.Date("2021-06-01"), as.Date("2025-02-17"))) +
  labs(x="Time (years)", y="Reassortment rate") +
  theme_minimal() +
  theme(legend.position="bottom") +
  scale_x_date(date_labels="%Y", date_breaks="1 year") +
  facet_wrap(name~., ncol=3, scales="free_y") +
  scale_alpha(guide=F) +
  ggtitle("Reassortment rate over time for HPAI and LPAI") +
  theme(plot.title = element_text(hjust = 0.5))
