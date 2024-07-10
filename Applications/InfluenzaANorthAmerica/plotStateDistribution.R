library(stringr)
library(seqinr)
library(ggplot2)
library("colorblindr")
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

h1n1color="#d95f02"
h3n2color="#1b9e77"

states <- c(
  "Alabama" = 5024279,
  "Arizona" = 7151502,
  "Arkansas" = 3011524,
  "California" = 39538223,
  "Colorado" = 5773714,
  "Connecticut" = 3605944,
  "Delaware" = 989948,
  "Florida" = 21538187,
  "Georgia" = 10711908,
  "Idaho" = 1839106,
  "Illinois" = 12812508,
  "Indiana" = 6785528,
  "Iowa" = 3190369,
  "Kansas" = 2937880,
  "Kentucky" = 4505836,
  "Louisiana" = 4657757,
  "Maine" = 1362359,
  "Maryland" = 6177224,
  "Massachusetts" = 7029917,
  "Michigan" = 10077331,
  "Minnesota" = 5706494,
  "Mississippi" = 2961279,
  "Missouri" = 6154913,
  "Montana" = 1084225,
  "Nebraska" = 1961504,
  "Nevada" = 3104614,
  "New_Hampshire" = 1377529,
  "New_Jersey" = 9288994,
  "New_Mexico" = 2117522,
  "New_York" = 20201249,
  "North_Carolina" = 10439388,
  "North_Dakota" = 779094,
  "Ohio" = 11799448,
  "Oklahoma" = 3959353,
  "Oregon" = 4237256,
  "Pennsylvania" = 13002700,
  "Rhode_Island" = 1097379,
  "South_Carolina" = 5118425,
  "South_Dakota" = 886667,
  "Tennessee" = 6916897,
  "Texas" = 29145505,
  "Utah" = 3271616,
  "Vermont" = 643077,
  "Virginia" = 8631393,
  "Washington" = 7693612,
  "West_Virginia" = 1793716,
  "Wisconsin" = 5893718,
  "Wyoming" = 576851,
  "District_Of_Columbia" = 712816
)

# read in the state ILI WHO_NREVSS_Clinical_Labs.csv in /cdcData/state, skip the first
state_ili = read.csv("./cdcData/states/WHO_NREVSS_Public_Health_Labs.csv", header=TRUE, sep=",", skip=1)

state_ili$Year = gsub("Season ", "", state_ili$SEASON_DESCRIPTION)
for (i in 1:nrow(state_ili)){
  state_ili$Year[i] = str_split(state_ili$Year[i], "-")[[1]][1]
}

# get all *HA.fasta files in the xml folder and read them in to count to number of samples per state
fasta_files = list.files("xmls", pattern="*HA.fasta", full.names=TRUE)
sample_counts = data.frame()
ratio = data.frame()
total = data.frame()
case_ratio = data.frame()

for (file in fasta_files) {
  # initialize a vector to keep track of the number of samples per state
  samples = c()
  for (state in names(states)) {
    samples[state] = 0
  }
  # Read in the fasta file
  seqs = read.fasta(file)
  # Get the US state for each sample from the header
  for (i in 1:length(seqs)) {
    # Get the header name
    header = names(seqs)[i]
    # split on / and get the second group, as the state
    state = str_split(header, "\\/")[[1]][2]
    state = gsub("_City", "", state)
    if (grepl("Human", state)){
      next
    }
    # Add one to the count for that state
    samples[state] = samples[state] + 1
  }
  # get the virus and year from the file name splitting on _
  tmp = str_split(file, "\\/")[[1]][2]
  tmp2 = str_split(tmp, "_")[[1]]
  # keep track of them as a dataframe using the state name as header
  sample_counts = rbind(sample_counts, data.frame( virus = tmp2[1], year = tmp2[2],  t(samples)))
  ratio_vals =c()
  totcount = c()
  cc = c()
  for (state in names(states)) {
    ratio_vals[state] = 0
    totcount[state] = 0
    cc[state] = 0
  }
  # for each state and season, get the corresponding number of cases
  for (state in names(states)) {
    # get the number of cases for that state and season
    if (grepl("H1N1", tmp2[1])){
      cases = as.numeric(state_ili[state_ili$REGION == state & state_ili$Year == tmp2[2], "A..2009.H1N1."])
    }else{
      cases = as.numeric(state_ili[state_ili$REGION == gsub("_"," ",state) & state_ili$Year == tmp2[2], "A..H3."])
    }
    # check if numeric(0) or NA
    if (length(cases) == 0 || is.na(cases)){
      cases = 0
    }
    # divide by the size of the state
    ratio_vals[state] = cases/states[state]
    totcount[state] = cases
    cc[state] = samples[state]/(cases+samples[state])
    if (cc[state] == Inf || is.na(cc[state])){
      cc[state] = 0
    }
  }
  
  # calculate the ratio of cases to samples
  ratio = rbind(ratio, data.frame( virus = tmp2[1], year = tmp2[2], t(ratio_vals)))
  total = rbind(total, data.frame( virus = tmp2[1], year = tmp2[2], t(totcount)))
  case_ratio = rbind(case_ratio, data.frame( virus = tmp2[1], year = tmp2[2], t(cc)))
}

# for each virus and season, compute the weighted average prevalence by weighting
# the prevalence of each state by the number of samples from that state
prev_norm = data.frame()
for (virus in unique(ratio$virus)){
  for (year in unique(ratio$year)){
    # get the ratio for that virus and year
    tmp = ratio[ratio$virus == virus & ratio$year == year,]
    # get the number of samples for that virus and year
    tmp2 = sample_counts[sample_counts$virus == virus & sample_counts$year == year,]
    # get the weighted average
    vals = c()
    vals1 = c()
    for (i in seq(3, ncol(tmp))){
      vals1 = c(vals1, as.numeric(tmp[,i]))
      vals = c(vals,  as.numeric(tmp[,i]) * as.numeric(tmp2[,i]))
    }
    weighted_avg = sum(vals)/sum(as.numeric(tmp2[,3:ncol(tmp2)]))
    # compute the total number of cases
    tot = total[total$virus == virus & total$year == year,]
    tot = sum(as.numeric(tot[,3:ncol(tot)]))
    
    # get the standard deviation of the casses to samples ratio
    cc = vals1/sum(vals)
    sd = sqrt(sum((cc - mean(cc))^2)/length(cc))
    prev_norm = rbind(prev_norm, data.frame(virus = virus, year = year, weighted = weighted_avg, total = tot, sd=sd, mean=mean(cc)))
  }
}

# print to file
write.csv(prev_norm, file = "weighted_prevalence.csv", row.names = FALSE)

#plot total vs. weighted suing ggplot
p = ggplot(prev_norm, aes(x=total, y=weighted)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x -1, color='black') +
  geom_point(aes(color=virus)) +
  geom_text(aes(label = year), hjust = 0, vjust = 0) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Total number of cases", y = "Weighted average prevalence") +
  ggtitle("Total number of cases vs. weighted average prevalence") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(h1n1color, h3n2color)) 
plot(p)
ggsave("./../../Figures/weighted_prevalence.pdf", width = 6, height = 6)


# for each season, compute the sd and mean using the data in prev_norm


# add an offset of +1 for H1N1 and -1 for H3N2 to case_ratio
prev_norm$offset = -1
prev_norm$offset[case_ratio$virus == "H1N1"] = 1

p = ggplot(prev_norm, aes(x=as.numeric(year)+offset*0.1, y=mean, color=virus)) +
  geom_point(aes(color=virus)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Standard deviation") +
  ggtitle("Standard deviation of cases to samples ratio") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(-0.01, 0.01)) +
  scale_color_manual(values = c(h1n1color, h3n2color))
plot(p)

