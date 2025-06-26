
# ==== 1) User settings: point these at your files ====
log_file   <- "~/Downloads/LPAI_PB2.log"
trees_file <- "~/Downloads/LPAI_PB2-lpai_pb2.trees"

library(ape)       # read.nexus(), branching.times()
library(ggplot2)   # plotting
library(dplyr)     # data handling
library(phytools)

log_df <- read.table(log_file, header=TRUE, sep="\t", comment.char="#")

group_cols <- grep("^GroupSizes\\.", names(log_df), value=TRUE)
pop_cols   <- grep("^PopSizes\\.",   names(log_df), value=TRUE)
stopifnot(length(group_cols)>0, length(pop_cols)>0)

# sort by their numeric suffix
group_idx  <- as.numeric(sub("^GroupSizes\\.", "", group_cols))
group_cols <- group_cols[order(group_idx)]
pop_idx    <- as.numeric(sub("^PopSizes\\.",   "", pop_cols))
pop_cols   <- pop_cols[order(pop_idx)]


# 4) COPY .trees → tmp, then append END; to that tmp only
tmp_trees <- tempfile(pattern="bicep_trees_", fileext=".trees")
file.copy(trees_file, tmp_trees, overwrite = TRUE)

# append END; to the **copy** only
cat("\nEND;\n", file = tmp_trees, append = TRUE)

# now safe to read
treelist <- read.nexus(tmp_trees)
tree     <- treelist[[1]]

# ==== 5) get sampling times from tip labels ====
# tip labels must be like "ID|foo|2025.123"
parts    <- strsplit(tree$tip.label, "\\|")
tip_times <- sapply(parts, function(x){
  if(length(x)<3) stop("Tip label '", paste(x,collapse="|"),
                       "' doesn’t have a 3rd '|' field")
  as.Date(x[3])
})
if(any(is.na(tip_times)))
  stop("Non‐numeric or missing times in some tip labels.")

max_time         <- max(tip_times, na.rm=TRUE)
min_time          <- min(tip_times, na.rm=TRUE)
time_diff = (max_time - min_time)/365

data = data.frame()
# use 100 samples to iterate over between 100 and length(treelsit)

for (i in floor(seq(length(treelist)/10, length(treelist), length.out=20))) {
  tree <- treelist[[i]]
  
  # 1) get the node‐heights *above the root* for every edge
  nh <- nodeHeights(tree)   # nedge×2: [parentHeight, childHeight]
  
  # 2) pick only the edges whose *child* is an internal node
  Ntip <- length(tree$tip.label)
  internal_edges <- which(tree$edge[,2] > Ntip)
  
  # 3) extract the *childHeight* for those edges → exactly one height per coalescent event
  node_heights <- nh[internal_edges, 2]
  
  # 4) if you need the tree‐height (i.e. the maximum tip‐height) for shifting into “time before sampling”:
  tree_height <- max(nh[,2])   # this sees all childHeights—including tips

  # 3) time before sampling for each node
  coal_times    <- tree_height - node_heights
  coal_times = c(tree_height, coal_times) 
  
  # build sorted times and intervals:
  all_times <- sort(coal_times, decreasing=F)

  # posterior‐mean for each interval
  k_est       <- as.numeric(log_df[i, group_cols])
  PopSize_est <- as.numeric(log_df[i, pop_cols])
  
  df = data.frame()
  curr_time = 0
  curr_index = 0
  for (j in seq(1, length(k_est))){
    next_time = all_times[curr_index+k_est[j]]
    df <- rbind(df, data.frame(
      start_time = curr_time,
      end_time   = next_time,
      PopSize    = PopSize_est[j]
    ))
    curr_time = next_time
    curr_index = curr_index + k_est[j]
  }
  
  
  # keep track of the Pop size at these time points in data
  for (time in seq(0,time_diff,0.1)) {
    # find the index for which time is between start_time and end_time
    index = which(df$start_time <= time & df$end_time > time)    

    data <- rbind(data, data.frame(
      time = time,
      PopSize = df[index, "PopSize"]
    ))
  }
}
# calculate the 95 % hpod for each time point in data
data_plot <- data.frame()
for (time in seq(0,time_diff,0.1)) {
  index = which(data$time == time)
  if (length(index) > 0) {
    upper = quantile(data$PopSize[index], 0.975)
    lower = quantile(data$PopSize[index], 0.025)
    mean = mean(data$PopSize[index])
    data_plot <- rbind(data_plot, data.frame(
      time = time,
      upper = upper,
      lower = lower,
      mean = mean
    ))
  }
}
data_plot$time = as.Date(max_time-data_plot$time*365)

# save the dataplot to a file in ./tables/
write.csv(data_plot, "./tables/LPAI_PB2_population_size.csv", row.names=FALSE)



# read in positive cases
cases = read.csv("./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv")
cases$date = as.Date(cases$Date_Collected, format="%Y-%m-%d")

# loop over all days
min_date = min(cases$date)
max_date = max(cases$date)
smoothed_case_data = data.frame()
for (d in seq(min_date+23, max_date-23, by="day")){
  # get all instances within that time window
  window = cases[cases$date >= d-23 & cases$date <= d+23, ]
  # get how many instances of Final_IAV are Detected
  total_AIV = sum(window$Final_IAV == "Detected", na.rm=TRUE)
  # get the number of high path cases 
  total_HPAI = sum(window$Final_H5 == "Detected" &  window$Final_Pathogenicity == "High Path AI", na.rm=TRUE)
  smoothed_case_data = rbind(smoothed_case_data, data.frame(
    date = d,
    positivity = (total_AIV-total_HPAI)/nrow(window),
    type = "LPAI"
  ))
  # smoothed_case_data = rbind(smoothed_case_data, data.frame(
  #   date = d,
  #   positivity = total_HPAI/nrow(window),
  #   type = "HPAI"
  # ))
}
smoothed_case_data$date = as.Date(smoothed_case_data$date, format="%Y-%m-%d")


p = ggplot(data_plot, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line() +
  labs(x="Time before sampling (years)", y="Effective Population Size") +
  theme_minimal() +
  ggtitle("LPAI PB2 Population Size Over Time") +
  geom_line(data=smoothed_case_data, aes(x=date, y=positivity*500, color=type, group=type), method=NA, size=0.5) +
  # add second axis
  scale_y_continuous(sec.axis = sec_axis(~./500, name = "LPAI test positivity")) +
  scale_x_date()+
  # scale_y_log10() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p, filename = "../../Figures/LPAI_PB2_population_size.pdf", width=8, height=4)

