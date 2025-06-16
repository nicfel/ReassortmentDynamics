library(stringr)
library(seqinr)
library(ggplot2)
library(treeio)
library(ggtree)

# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in mcc.tree using treeio beast.tree
mcc.tree <- read.beast("./mcc.tree")

mcc.tree@data$between = as.numeric(as.numeric(mcc.tree@data$`HPAI+LPAI`)>0.5)

mcc.tree@data$between = as.numeric(as.numeric(mcc.tree@data$`HPAI+LPAI`)>0.5)*2
mcc.tree@data$between[is.na(mcc.tree@data$between)] = 0

# plot using ggtree, color by HPAI+LPAI
p = ggtree(mcc.tree, aes(color=between), mrsd=as.Date("2024-04-29")) +
  theme_tree2()+
  geom_tippoint(aes(fill=type), size=1.25, shape=21) +
  coord_cartesian(xlim=c(2021, 2025)) +
  scale_color_continuous(limits=c(0, 1),  low="grey40", high="#ff7f00", name="reassortment events\nbetween\nHigh & Low Path.") +
  scale_fill_manual(values=c("#4daf4a", "#377eb8"), name="") 

plot(p)
ggsave("HPAILPAItree.pdf", p, width=8, height=6, units="in", dpi=300)