library(ggplot2)
library(ggpubr)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# set colors of structured with and without to ['#1b9e77','#d95f02','#7570b3']
sims_colors = c("without"="#1b9e77", "with"="#d95f02", "structured"="#7570b3")


t1 = read.table("./waitTimes.tsv", header = FALSE, sep = "\t")
t2 = read.table("./structWaitTimes.tsv", header = FALSE, sep = "\t")
t2$V1 = t2$V1 + max(t1$V1)

t = rbind(t1, t2)

t$superspreading = "without"
t[t$V1>500, "superspreading"] = "with"
t[t$V1>1000, "superspreading"] = "structured"
t$V3 = as.numeric(t$V2)

# simulate exponential random numbers with mean 1 for comparison
n = 100000
x = rexp(n, rate = 1)
ks_data = data.frame()
for (approach in unique(t$superspreading)) {
  subset = t[t$superspreading == approach, "V2"]
  # get the KS distance of subset from the exponential distribution for increasing numbers of samples
  vals = round(seq(n/10,length(subset),length.out = 200))
  ks = sapply(vals, function(i) ks.test(subset[1:i], x)$statistic)
  # make a dataframe wit the ks values
  ks_data = rbind(ks_data, data.frame(n = vals, ks = ks, superspreading = approach))
}

# make 4 plots, first, make a histogram of distributions for with super spreading to 
# the exponential distribution, scale the y axis to be relative frequency
expo_df = data.frame(x = x)

p1 = ggplot(t[t$superspreading == "without",], aes(x = V2)) + 
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), binwidth=0.1, geom="step", position="identity") +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), binwidth = 0.1, geom="step", position="identity") +
  theme_minimal()+
  # geom_density(aes(y = ..count..), fill = "blue", alpha = 0.5) + 
  ggtitle("Without Super Spreading") + 
  # put legend to the top right within the plot
  theme(legend.position = c(0.7, 0.7))+
  coord_cartesian(xlim = c(0, 3))+
  xlab("Wait Time") + 
  ylab("Frequency")
plot(p1)

p2 = ggplot(t[t$superspreading == "with",], aes(x = V2)) + 
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), binwidth=0.1, geom="step", position="identity") +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), binwidth = 0.1, geom="step", position="identity") +
  theme_minimal()+
  # geom_density(aes(y = ..count..), fill = "blue", alpha = 0.5) + 
  ggtitle("With Super Spreading") + 
  scale_color_viridis_d(name="")+
  coord_cartesian(xlim = c(0, 3))+
  xlab("Wait Time") + 
  ylab("Frequency")
plot(p2)

p3 = ggplot(t[t$superspreading == "structured",], aes(x = V2)) + 
  # geom_histogram(data=expo_df, aes(x=x,y = ..density.., fill="exponential distribution"), binwidth = 0.1, alpha = 0.25) +
  # geom_histogram(aes(y = ..density..,  fill = "wait time distribution"), binwidth = 0.1, alpha = 0.25) +
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), binwidth=0.1, geom="step", position="identity") +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), binwidth = 0.1, geom="step", position="identity") +
  theme_minimal()+
  ggtitle("With 20 populations and superspreading") + 
  scale_color_viridis_d(name="")+
  xlab("Wait Time") + 
  ylab("Frequency")
plot(p3)


# plot the ks distance
p4 = ggplot(ks_data, aes(x = n, y = ks, color=superspreading)) + 
  geom_line() + 
  theme_minimal() + 
  ggtitle("KS distance from exponential distribution") + 
  xlab("Number of samples") + 
  ylab("KS distance")+
  # remove x axix
  scale_color_manual(values = sims_colors)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  # put legend to the top right within the plot
  theme(legend.position = c(0.9, 0.9))+
  scale_x_log10()
plot(p4)



# plot them all together using ggarrange
p = ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
plot(p)
ggsave("./../../Figures/WaitTimeDistribution.pdf", p, width = 9, height = 6)



