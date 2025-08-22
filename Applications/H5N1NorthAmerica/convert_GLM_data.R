library(stringr)
library(seqinr)
library(ggplot2)
library("colorblindr")
library(ggtree)
library(treeio)
library(ggnewscale)
library(cowplot)

rm(list=ls())

clades = c("B3.13", "D1.1")
rate_shift_str = '0 0.105936073059361 0.211872146118721 0.317808219178082 0.423744292237443 0.529680365296804 0.635616438356164 0.741552511415525 0.847488584474886 0.953424657534247 1.05936073059361 1.16529680365297 1.27123287671233 1.37716894977169 1.48310502283105 1.58904109589041 1.69497716894977 1.80091324200913 1.90684931506849 2.01278538812785 2.11872146118721 2.22465753424658 2.33059360730594 2.4365296803653 2.54246575342466 2.64840182648402 2.75433789954338 2.86027397260274 2.9662100456621 3.07214611872146 3.17808219178082 3.28401826484018 3.38995433789954 3.4958904109589 3.60182648401827 3.70776255707763 3.81369863013699 3.91963470319635 4.02557077625571 4.13150684931507'
rate_shifts = as.numeric(strsplit(rate_shift_str, " ")[[1]])

system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ","-burnin 20 -log ./out/HPAI_HLHxNx.glm.rep*.trees -o ./combined/HPAI_HLHxNx.glm.trees"))
system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher ReassortmentNetworkSummarize -burnin 0 -followSegment 0  -positions MCC  ./combined/HPAI_HLHxNx.glm.trees ./combined/HPAI_HLHxNx.glm.tree"))
system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner -burnin 20 -log ./out/HPAI_HLHxNx.glm.rep*.log -o ./combined/HPAI_HLHxNx.glm.log"))

system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher GetCladeHeightsFromNetwork",
             "-burnin 0 -tree ./combined/HPAI_HLHxNx.glm.trees -clade ./tables/cow_clade.csv -out ./combined/b113_glm.tsv"))
system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher GetCladeHeightsFromNetwork",
             "-burnin 0 -tree ./combined/HPAI_HLHxNx.glm.trees -clade ./tables/d11.csv -out ./combined/d11_glm.tsv"))

segment_order = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

for (s in seq(1,length(segment_order), 1)){
  # get the segment name
  segment = segment_order[s]
  # run the applauncher to mark the clades for this segment
  system(paste0("/Applications/BEAST\\ 2.7.7/bin/applauncher MarkCladesFromCladeFile ",
               "-burnin 0 -followSegment ", s-1, " -printSegment ", s-1,
               " -tree ./combined/HPAI_HLHxNx.glm.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HPAI_HLHxNx.glm.", segment, ".trees"))
}
# das

# Read in the clade heights
clade_cow_heights <- read.csv("./combined/b113_glm.tsv", sep="\t")
clade_d11_heights <- read.csv("./combined/d11_glm.tsv", sep="\t")
log_file <- read.csv("./combined/HPAI_HLHxNx.glm.log", sep="\t")
data = data.frame();


for (cl in clades){
  if (cl == "B3.13"){
    clade_heights <- clade_cow_heights
  }else if (cl == "D1.1"){
    clade_heights <- clade_d11_heights
  }
  
  # define the rate shifts values
  #loop over the posterior
  no_event_probs = c()
  for (l in seq(1, min(nrow(clade_heights), nrow(log_file)), 1)) {
    # get the timings of the HA segment
    if (cl == "B3.13"){
      min_time = clade_heights[l, 2]
      max_time = min_time+0.5
    }else if (cl == "D1.1"){
      min_time = clade_heights[l, 2]
      max_time = min_time+0.5
    }
    
    
    first_interval = which(rate_shifts <= min_time)[length(which(rate_shifts <= min_time))]
    last_interval = which(rate_shifts <= max_time)[length(which(rate_shifts <= max_time))]
    
    curr_time = min_time
    weighted = 0.0

    for (i in seq(first_interval, last_interval)) {
      if (i > length(rate_shifts)) {
        stop("rate shifts out of bounds")
      }
      next_time = min(rate_shifts[i + 1], max_time)
      
      # get the reassortment rates of this interval at the beginning and end
      r_start = log_file[l, paste0("reassortmentRate.", i)]
      r_end = log_file[l, paste0("reassortmentRate.", i + 1)]

      # calculate the growth rate for this interval  
      growth = (r_start - r_end)/(rate_shifts[i+1]-rate_shifts[i]);
      
      timediff1 = curr_time - rate_shifts[i]
      timediff2 = next_time - rate_shifts[i]
      
      if (growth == 0.0) {
        weighted = weighted +  (next_time - curr_time) * Math.exp(rates[i]);
      } else {
        weighted <- weighted + exp(r_start)/(-growth) * (
          exp(-growth * timediff2) - exp(-growth * timediff1)
        )
      }
      curr_time = next_time;
    }
    
    # Compute probability of no event over interval
    data = rbind(data, data.frame(
      no_event_prob = 1-(1-0.5^(8-1))*exp(-weighted),
      min_time = min_time,
      max_time = max_time,
      mean_rate = weighted / (max_time - min_time),
      weighted = weighted,
      clade = cl
    ))
  }
}

# print data to csv file
write.csv(data, "./combined/HPAI_HLHxNx.glm.cladeprobs.csv", row.names = FALSE)

