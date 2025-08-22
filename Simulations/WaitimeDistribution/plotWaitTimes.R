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

# Define consistent color scheme
wait_time_color <- "#d95f02"  # Orange for wait time distribution
exponential_color <- "#1b9e77"  # Teal for exponential distribution

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
expo_df = data.frame(x = x)

p1 = ggplot(t[t$superspreading == "without",], aes(x = V2)) + 
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), 
           binwidth=0.1, geom="step", position="identity", size=0.8) +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), 
           binwidth = 0.1, geom="step", position="identity", size=0.8) +
  scale_color_manual(values = c("exponential distribution" = exponential_color,
                                "wait time distribution" = wait_time_color),
                     name = "") +
  theme_minimal() +
  ggtitle("SIR") + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(0, 3)) +
  xlab("Wait Time") + 
  ylab("Density")

p2 = ggplot(t[t$superspreading == "with",], aes(x = V2)) + 
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), 
           binwidth=0.1, geom="step", position="identity", size=0.8) +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), 
           binwidth = 0.1, geom="step", position="identity", size=0.8) +
  scale_color_manual(values = c("exponential distribution" = exponential_color,
                                "wait time distribution" = wait_time_color),
                     name = "") +
  theme_minimal() +
  ggtitle("SIR with superspreading") + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(0, 3)) +
  xlab("Wait Time") + 
  ylab("Density")

p3 = ggplot(t[t$superspreading == "structured",], aes(x = V2)) + 
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), 
           binwidth=0.1, geom="step", position="identity", size=0.8) +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), 
           binwidth = 0.1, geom="step", position="identity", size=0.8) +
  scale_color_manual(values = c("exponential distribution" = exponential_color,
                                "wait time distribution" = wait_time_color),
                     name = "") +
  theme_minimal() +
  ggtitle("20 state SIR with superspreading") + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("Wait Time") + 
  ylab("Density")

# Extract legend from p1 (before removing its legend)
temp_plot_for_legend <- ggplot(t[t$superspreading == "without",], aes(x = V2)) + 
  stat_bin(data=expo_df, aes(x=x, y=..density.., color="exponential distribution"), 
           binwidth=0.1, geom="step", position="identity", size=0.8) +
  stat_bin(aes(y = ..density.., color = "wait time distribution"), 
           binwidth = 0.1, geom="step", position="identity", size=0.8) +
  scale_color_manual(values = c("exponential distribution" = exponential_color,
                                "wait time distribution" = wait_time_color),
                     name = "") +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"))

library("cowplot")
# Extract just the legend
legend_only <- get_legend(temp_plot_for_legend)

p = plot_grid(p1, p2, p3, legend_only, labels = c("A", "B", "C", ""), ncol = 2)

plot(p)
ggsave("./../../Figures/WaitTimeDistribution.pdf", p, width = 6, height = 4)