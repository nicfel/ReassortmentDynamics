library(stringr)
library(seqinr)
library(ggplot2)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# print the working directory
getwd()
log.files <- list.files("./out/", pattern="*.log", full.names=TRUE)


# read in and combine the three constant and the three varibale files
# read in the constant files
constant1 = read.table("./out//H5N1_wgs.constant.rep0.log",header=TRUE)
constant2 = read.table("./out//H5N1_wgs.constant.rep1.log",header=TRUE)
constant3 = read.table("./out//H5N1_wgs.constant.rep2.log",header=TRUE)
# combine the files after a 10 % burnin
constant = rbind(constant1[-c(1:round(0.1*nrow(constant1))),],
                 constant2[-c(1:round(0.1*nrow(constant2))),],
                 constant3[-c(1:round(0.1*nrow(constant3))),])

# read in the variable files
variable1 = read.table("./out//H5N1_wgs.varying.rep0.log",header=TRUE)
variable2 = read.table("./out//H5N1_wgs.varying.rep1.log",header=TRUE)
variable3 = read.table("./out//H5N1_wgs.varying.rep2.log",header=TRUE)
# combine the files after a 10 % burnin
variable = rbind(variable1[-c(1:round(0.1*nrow(variable1))),],
                 variable2[-c(1:round(0.1*nrow(variable2))),],
                 variable3[-c(1:round(0.1*nrow(variable3))),])

# plot the reassortment rates 1 trhough 3 for variable and 1 for constant
data = data.frame(
  from=0, to=3,
  rate = median(constant$reassortmentRate.1),
  lower = quantile(constant$reassortmentRate.1, 0.025),
  upper = quantile(constant$reassortmentRate.1, 0.975),
  method="constant"
)

data = rbind(data,data.frame(
  from=0, to=1,
  rate = median(variable$reassortmentRate.1),
  lower = quantile(variable$reassortmentRate.1, 0.025),
  upper = quantile(variable$reassortmentRate.1, 0.975),
  method="variable"
))
data = rbind(data, data.frame(
  from=1, to=2,
  rate = median(variable$reassortmentRate.2),
  lower = quantile(variable$reassortmentRate.2, 0.025),
  upper = quantile(variable$reassortmentRate.2, 0.975),
  method="variable"
))
data = rbind(data, data.frame(
  from=2, to=3,
  rate = median(variable$reassortmentRate.3),
  lower = quantile(variable$reassortmentRate.3, 0.025),
  upper = quantile(variable$reassortmentRate.3, 0.975),
  method="variable"
))


data$from = as.Date("2024-02-28")-0.806392694063927*data$from*365
data$to = as.Date("2024-02-28")-0.806392694063927*data$to*365

p = ggplot(data, aes(x=from, xend=to, y=rate, yend=rate, color=method)) +
  geom_segment(size=1) +
  # add the upper and lower bounds and fill in the area between as segments
  geom_segment(aes(x=from, xend=from, y=lower, yend=lower), size=0.5) +
  geom_segment(aes(x=to, xend=to, y=upper, yend=upper), size=0.5) +
  # fill in the square between upper and lower
  geom_rect(aes(xmin=from, xmax=to, ymin=lower, ymax=upper, fill=method), alpha=0.2) +
  theme_minimal()
plot(p)




