library(stringr)
library(seqinr)
library(ggplot2)
library(cowplot)


# 1st part of plotDynamics.R
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)


segments = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

# Set the directory to the directory of the file
setwd("~/all_HPAI/")

# get all tree files with rep0 in out
tree.files = list.files("variable", pattern="rep0.trees", full.names=T)
# combine the tree file with the rep1 and rep2 using log combiner
# this is commented out because it takes a long time, but has to be done once
# TODO uncomment
for (i in 1:length(tree.files)) {
  system(
    paste("/Applications/BEAST\\ 2.7.6/bin/logcombiner -burnin 80 -log", tree.files[i],
          gsub("rep0", "rep1", tree.files[i]), gsub("rep0", "rep2", tree.files[i]),
          "-o",
          gsub(".rep0.trees", ".combined.trees", tree.files[i])))
}
# get all the combined trees
tree.files = list.files("variable", pattern="combined.trees", full.names=T)
# run them through the  app to create the MCC network
# TODO uncomment
for (i in 1:length(tree.files)) {
  system(
    paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentNetworkSummarizer",
          "-burnin 0 -positions none",
          "-followSegment 0",tree.files[i],
          gsub(".trees", ".tree", tree.files[i])))
  for (a in seq(0,7)){
    for (b in seq(a+1,7)){
      remove_segments = seq(0,7)
      # remove the entries ==a and ==b
      remove_segments = remove_segments[remove_segments!=a]
      remove_segments = remove_segments[remove_segments!=b]
      system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
             paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentNetworkSummarizer",
                   "-burnin 0 -positions none -removeSegments",paste(remove_segments, collapse=","),
                   "-followSegment", a, tree.files[i],
                   gsub(".trees", paste(".", a,"_",b, ".tree", sep=""), tree.files[i])))
      
      system(intern=T, show.output.on.console=F,ignore.stdout = T, ignore.stderr = T,
             paste("/Applications/BEAST\\ 2.7.6/bin/applauncher ReassortmentOverTime",
                   "-burnin 0 -removeSegments",paste(remove_segments, collapse=",") , tree.files[i],
                   gsub(".trees", paste(".", a,"_",b, ".reassortment.txt", sep=""), tree.files[i])))
      
    }
  }
}








