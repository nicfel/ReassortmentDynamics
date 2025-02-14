library(stringr)
library(seqinr)
library(ggplot2)
library(cowplot)

# 3rd part of plotDynamics.R

segments = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

# Set the directory to the directory of the file
setwd("~/all_HPAI/")

# get all the combined trees
tree.files = list.files("variable", pattern="combined.trees", full.names=T)
tree.files = list.files("variable", pattern=".tree", full.names=T)

# read in the reassortment.txt files, the first row is the run number, the second
# is the time, the third is 1/rate
reassortment = data.frame()
maxtime = 3
# for (i in 1:length(tree.files)) {
#   t = read.table(gsub(".tree", ".reassortment.txt", tree.files[i]), header=F, sep="\t")
#   for (run in unique(t$V1)) {
#     tmp = t[t$V1==run,]
#     timeFrom = c(0, tmp$V2)
#     name = strsplit(tree.files[i], split="\\.")[[1]]
#     reassortment = rbind(reassortment, data.frame(run=run,
#                                                   from = timeFrom[1:(length(timeFrom)-1)],
#                                                   to = tmp$V2,
#                                                   rate=1/tmp$V3, 
#                                                   file=name[[2]]))
#   }
# }

reassortmentfiles = list.files("variable/", pattern="reassortment.txt", full.names=T)
# remove the files with 7_7 and 7_8
reassortmentfiles = reassortmentfiles[!grepl("7_7", reassortmentfiles)]
reassortmentfiles = reassortmentfiles[!grepl("7_8", reassortmentfiles)]
counts = data.frame()
for (i in 1:length(reassortmentfiles)) {
  name = strsplit(reassortmentfiles[i], split="\\.")[[1]]
  if (length(name)==6){
    t = read.table(reassortmentfiles[i], header=F, sep="\t")
    fromto = strsplit(name[[4]], split="_")[[1]]
    # for each unique V0, compute how many instances of V2 <maxtime there are
    for (run in unique(t$V1)) {
      tmp = t[t$V1==run,]
      counts = rbind(counts, data.frame(run=run, 
                                        from=fromto[1], to=fromto[2], 
                                        count=sum(tmp$V2 < maxtime),
                                        file=name[[2]]))
      counts = rbind(counts, data.frame(run=run, 
                                        from=fromto[2], to=fromto[1], 
                                        count=sum(tmp$V2 < maxtime),
                                        file=name[[2]]))
      
    }
  }
}

segments = c("HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA")

colors = c('#fdbf6f', '#e31a1c', '#ff7f00', '#a6cee3', '#1f78b4', '#fb9a99', '#b2df8a', '#33a02c')

colors = c('HA'='#fdbf6f', 'NA'='#ff7f00','MP'='#a6cee3', 'NP'='#1f78b4', 'PA'='#66c2a4','PB1'='#2ca25f','PB2'='#006d2c', 'NS'='#fb9a99')
# plot a histogram for file=="constant" for the count for each combination of from to
#rename from and to from "0" "1".. "7" to segments
counts$froms = segments[as.numeric(counts$from)+1]
counts$tos = segments[as.numeric(counts$to)+1]
counts$froms = factor(counts$froms, levels=labels(colors))
counts$tos = factor(counts$tos, levels=labels(colors))

p = ggplot(counts, aes(x=count, fill=tos, color=tos)) + 
  geom_density(alpha=0.4) + 
  xlab("Reassortment events") + 
  facet_grid(froms~file ) +
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  ggtitle("Reassortment counts - HPAI+LPAI") + 
  theme_minimal() +
  theme(text= element_text(size = 18)) 
  
plot(p)
ggsave("Figures/H5N1_CoReassortmentDistr_HPAI_all.pdf", p, width = 11, height = 8.5, units = "in")


p <- ggplot(counts[counts$froms == "HA" & counts$file == "variable",], aes(x = count, fill = tos, color = tos)) + 
  geom_density(alpha = 0.4) + 
  xlab("Reassortment events") +
  scale_color_manual(name="", values = colors) +
  scale_fill_manual(name="", values = colors) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 40)) +
  theme(
    text = element_text(size = 18),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(), 
    axis.line.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom"
  ) + ggtitle("Reassortments with HA (HPAI+LPAI)")
  
  
plot(p)
ggsave("Figures/H5N1_HA_variable _CoReassortmentDistr_HPAI_all.pdf", p, width = 11, height = 8.5, units = "in")

# plot the mean counts as a tile plot
library(dplyr)
mean_counts <- counts %>%
  group_by(from, to, file) %>%
  summarise(mean_count = mean(count, na.rm = TRUE)) %>%
  ungroup()

mean_counts$froms = segments[as.numeric(mean_counts$from)+1]
mean_counts$tos = segments[as.numeric(mean_counts$to)+1]
mean_counts$froms = factor(mean_counts$froms, levels=labels(colors))
mean_counts$tos = factor(mean_counts$tos, levels=rev(labels(colors)))

p = ggplot(mean_counts, aes(x=froms, y=tos, fill=log(mean_count))) + 
  geom_tile() + facet_wrap(~file) + 
  scale_fill_viridis_c(name="log number\nevents")+
  xlab("") + 
  ylab("") + 
  ggtitle("CoReassortment mean number - HPAI+LPAI") + 
  theme_minimal() + 
  theme(text = element_text(size = 18))

plot(p)
ggsave("Figures/H5N1_CoReassortmentMean_HPAI_all.pdf", p,  width = 11, height = 8.5, units = "in")

# for each file in reassormtne compute a moving average from time =0 to time = 3
# read in the reassortment.txt files, the first row is the run number, the second
# is the time, the third is 1/rate
reassortment = data.frame()
maxtime = 3
 for (i in 1:length(tree.files)) {
   t = read.table(gsub(".tree", ".reassortment.txt", tree.files[i]), header=F, sep="\t")
   for (run in unique(t$V1)) {
     tmp = t[t$V1==run,]
     timeFrom = c(0, tmp$V2)
     name = strsplit(tree.files[i], split="\\.")[[1]]
     reassortment = rbind(reassortment, data.frame(run=run,
                                                   from = timeFrom[1:(length(timeFrom)-1)],
                                                   to = tmp$V2,
                                                   rate=1/tmp$V3, 
                                                   file=name[[2]]))
   }
 }


timepoints = seq(0, 3, 0.1)
rea_avg=data.frame()

# mrsi date can be found in script that calls XML (script for variable rates estimation)
# note the following function takes a long time 
for (f in unique(reassortment$file)) {
  tmp = reassortment[reassortment$file==f,]
  for (time in timepoints) {
    rate = tmp[tmp$from <= time & tmp$to >= time, "rate"]
    rea_avg = rbind(rea_avg, data.frame(from=time, to=time+0.1, fromtime=mrsi-j*max(time)/1000*365, rate=mean(rate), 
                                        lower=quantile(rate, 0.025), upper=quantile(rate, 0.975),
                                        file=f))
  }
}

# plot everything as a step function in ggplot
p = ggplot(rea_avg, aes(x=from, xend=to, y=rate, yend=rate)) + 
  geom_segment() +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) +
  facet_wrap(~file) +
  xlab("Time (years)") + 
  #coord_cartesian(ylim=c(0,4)) +
  ylab("Reassortment rate") + 
  ggtitle("Reassortment rate over time (HPAI_all)") + 
  theme_minimal() +
  theme(text = element_text(size = 18))

plot(p)
ggsave("Figures/H5N1_reassortment_rate_overtime_HPAI_all.pdf", p,  width = 11, height = 8.5, units = "in")




### export all the dataframes
counts
mean_counts

write.csv(rea_avg, "HPAI_all_variable_moving_avreage_reasstormentrate.csv")
