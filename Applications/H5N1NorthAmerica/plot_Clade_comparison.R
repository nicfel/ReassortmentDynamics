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
rate_shift_str = '0 1.25 2.5 3.75 5 6 9.5 13 16.5 20'
rate_shifts = as.numeric(strsplit(rate_shift_str, " ")[[1]])

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

rerun = F

for (isconstant in c(TRUE)){
  
  if (isconstant){
    if (rerun){
      # run log combined on the constant trees
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out2/HLHxNx.constant.rep*.trees -o ./combined/HLHxNx.constant.trees"))
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out2/HLHxNx.constant.rep*.trees -o ./combined/HLHxNx.constant.trees"))
      

      system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher ReassortmentNetworkSummarize -burnin 0 -followSegment 0  -positions MCC  ./combined/HLHxNx.constant.trees ./combined/HLHxNx.constant.tree"))
      # system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher ReassortmentNetworkSummarize -burnin 0 -followSegment 0  -positions MCC  ./combined/HPAI_HLHxNx.constant.trees ./combined/HPAI_HLHxNx.constant.tree"))
      
      # system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner -burnin 20 -log ./out2/HLHxNx.constant.rep*.log -o ./combined/HLHxNx.constant.log"))
      # system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner -burnin 20 -log ./out/LPAI_HLHxNx.constant.rep*.log -o ./combined/LPAI_HLHxNx.constant.log"))
      # system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner -burnin 20 -log ./out/HPAI_HLHxNx.constant.rep*.log -o ./combined/HPAI_HLHxNx.constant.log"))

      
      system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher MarkCladesFromCladeFile",
                   "-burnin 0 -followSegment 0 -tree ./combined/HLHxNx.constant.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.constant.clades.trees"))
      system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher MarkCladesFromCladeFile",
                   "-burnin 0 -followSegment 0 -printTable true -tree ./combined/HLHxNx.constant.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.constant.clades.tsv"))
      
      
      for (s in seq(1,length(segment_order), 1)){
        # get the segment name
        segment = segment_order[s]
        # run the applauncher to mark the clades for this segment
        system(paste0("/Applications/BEAST\\ 2.7.7/bin/applauncher MarkCladesFromCladeFile ",
                     "-burnin 0 -followSegment ", s-1, " -printSegment ", s-1,
                     " -tree ./combined/HLHxNx.constant.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.constant.", segment, ".trees"))
        system(paste0("/Applications/BEAST\\ 2.7.7/bin/treeannotator ",
                     "-burnin 0 -height keep ./combined/HLHxNx.constant.", segment, ".trees  ./combined/HLHxNx.constant.", segment, ".tree"))
      }
    }
    
    
    # read in the log file for HLHxNx
    log_file <- read.csv("./combined/HLHxNx.constant.log", sep="\t")
    # log_file_lpai <- read.csv("./combined/LPAI_HLHxNx.constant.log", sep="\t")
    # log_file_hpai <- read.csv("./combined/HPAI_HLHxNx.constant.log", sep="\t")
  }else{
    
    if (rerun){
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out/HLHxNx.dependent.rep*.trees -o ./combined/HLHxNx.dependent.trees"))
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out/HPAI_HLHxNx.dependent.rep*.trees -o ./combined/HPAI_HLHxNx.dependent.trees"))
      
      system(paste("/Applications/BEAST\\ 2.7.7/bin/applauncher ReassortmentNetworkSummarize",
                   "-burnin 0 -followSegment 0  -positions MCC  ./combined/HPAI_HLHxNx.dependent.trees ./combined/HPAI_HLHxNx.dependent.tree"))
      
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out/HLHxNx.dependent.rep*.log -o ./combined/HLHxNx.dependent.log"))
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out/LPAI_HLHxNx.dependent.rep*.log -o ./combined/LPAI_HLHxNx.dependent.log"))
      system(paste("/Applications/BEAST\\ 2.7.7/bin/logcombiner ",
                   "-burnin 20 -log ./out/HPAI_HLHxNx.dependent.rep*.log -o ./combined/HPAI_HLHxNx.dependent.log"))
      

      # for (s in seq(1,length(segment_order), 1)){
      #   # get the segment name
      #   segment = segment_order[s]
      #   # run the applauncher to mark the clades for this segment
      #   system(paste0("/Applications/BEAST\\ 2.7.7/bin/applauncher MarkCladesFromCladeFile ",
      #                "-burnin 0 -followSegment ", s-1, " -printSegment ", s-1,
      #                " -tree ./combined/HLHxNx.dependent.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.dependent.", segment, ".trees"))
      #   system(paste0("/Applications/BEAST\\ 2.7.7/bin/treeannotator ",
      #                "-burnin 0 -height keep ./combined/HLHxNx.dependent.", segment, ".trees  ./combined/HLHxNx.dependent.", segment, ".tree"))
      # }
    }

    # read in the log file for HLHxNx
    log_file <- read.csv("./combined/HLHxNx.dependent.log", sep="\t")
    # log_file_lpai <- read.csv("./combined/LPAI_HLHxNx.dependent.log", sep="\t")
    # log_file_hpai <- read.csv("./combined/HPAI_HLHxNx.dependent.log", sep="\t")
  }
        

  # # put your 3 data.frames into a named list
  # df_list <- list(
  #   both = log_file
  # )
  # 
  # for (df_name in names(df_list)) {
  #   lf <- df_list[[df_name]]
  #   # choose the appropriate mrsi vector
  #   mrsi_tmp <- if (df_name == "lpai") mrsi_lpai else mrsi
  #   for (i in seq(1, length(rate_shifts))){
  #     # get the reassortment rate at this time point
  #     rate = lf[, paste0("InfectedToRho.", i)]
  #     if (!isconstant){
  #       rate = rate+lf[, paste0("logNe.", i)]
  #     }
  # 
  #     # get all the quantile from 0.05 to 1.0
  #     for (q in seq(0.05, 1.0, 0.05)){
  #       upper = quantile(rate, 1-q/2)
  #       lower = quantile(rate, q/2)
  #       if (q==1){
  #         lower = lower+0.03
  #         upper = upper-0.03
  #       }
  #       reassortment_rate = rbind(reassortment_rate, data.frame(
  #         time = mrsi_tmp - rate_shifts[i]*365,
  #         quantile = q,
  #         upper = upper,
  #         lower = lower,
  #         isconstant = isconstant,
  #         name = df_name
  #       ))
  #     }
  #   }
  # }
}

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
  smoothed_case_data = rbind(smoothed_case_data, data.frame(
    date = d,
    positivity = total_HPAI/nrow(window),
    type = "HPAI"
  ))
}
smoothed_case_data$date = as.Date(smoothed_case_data$date, format="%Y-%m-%d")


# start new variable, constant rho
reassortment_rate$method = "constant rho"

# 1) Methods (3 distinct hues from Set1)
methods_colors <- c(
  "constant"= "#E41A1C",  # red
  "constant rho"= "#377EB8",  # blue
  "dependent rho(t)"= "#4DAF4A"   # green
)

reassortment_rate$alpha = 0.2
# for q==1 set it to 1
reassortment_rate$alpha[reassortment_rate$quantile == 1] = 1.0
# start new variable, constant rho

y_val = log(0.03)
y_off = 0.2

for (ind in c(TRUE)){
  # plot the reassortment rate over time
  smoothed_case_data$method=NA
  smoothed_case_data$quantile=NA
  lineage_colors <- c(
    hpai    = "#E41A1C",  # red
    lpai    = "#377EB8",  # blue
    both  = "#4DAF4A"  # green
  )
  smoothed_case_data$name =NA
  # p = ggplot(reassortment_rate[reassortment_rate$isconstant == ind, ], aes(x=time, y=upper, group=interaction(quantile, name), fill=name)) +
  #   geom_ribbon(aes(ymin=lower, ymax=upper, alpha=alpha)) +
  #   coord_cartesian(xlim=c(as.Date("2021-09-01"), mrsi), ylim=c(-4, log(5))) +
  #   scale_alpha(guide=F) +
  #   scale_fill_manual(values=lineage_colors, name="Dataset") +
  #   scale_color_manual(values=clade_colors, name="Median Emergence Time of Clade:") +
  #   new_scale_color()+
  #   # geom_line(data=smoothed_case_data, aes(x=date, y=log(positivity), color=type, group=type), method=NA, size=0.5) +
  #   scale_color_manual(values=c("HPAI"="#E41A1C", "LPAI"="#377EB8"), name="Type") +
  #   
  #   scale_y_continuous(breaks=c(log(0.05), log(0.1), log(0.2), log(0.4), log(0.8), log(1.6), log(3.2)), 
  #                      labels=c("0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "3.2")) +
  #   theme_minimal() + 
  #   theme(legend.position = "top")+
  #   xlab("Time") +
  #   ylab("Reassortment rate") 
  # plot(p)
  # if (ind){
  #   ggsave(p, filename="../../Figures/h5n1_reassortment_rate_constant.pdf", width=6, height=2.5)
  # }else{
  #   ggsave(p, filename="../../Figures/h5n1_reassortment_rate_dependent.pdf", width=6, height=2.5)
  # }



  
  # read in ./combined/HLHxNx.constant.HA.tree
  if (ind){
    tree = read.beast("./combined/HLHxNx.constant.HA.tree")
  }else{
    tree = read.beast("./combined/HLHxNx.dependent.HA.tree")
  }

  # plot the tree, coloring by the trait HPAI+LPAI
  tree@data$loc = as.numeric(tree@data$`HPAI+LPAI`)
  # replace all NA with 0
  tree@data$loc[is.na(tree@data$loc)] = 0
  tree@data$loc[tree@data$loc>1] = 1
  tree@data$iscow = "bird"
  # which tip‐indices contain "cow"?
  cow_tips <- which(grepl("cow", tree@phylo$tip.label, ignore.case = TRUE))
  # now flag those rows in @data whose node is in cow_tips
  tree@data$iscow[tree@data$node %in% cow_tips ] <- "cow"
  # convert mrsi to decimal year
  mrsi_dec <- as.numeric(mrsi - as.Date("1970-01-01")) / 365.25 + 1970
  
  # now plot
  p_tree <- ggtree(tree, aes(color = loc), mrsd = mrsi, size = 0.2) +
    theme_tree2() +
    # add cow tips as filled circles
    coord_cartesian(xlim=c(2020, mrsi_dec+0.05)) +
    # add node labels for posterior support
    scale_color_gradient2(low="black", mid="#F93946", high="#F93946", midpoint=0.7, 
                          name="posterior support for\nHPAI LPAI reassortment",
                          breaks=c(0,0.5,1), labels=c("0", "0.5", "1+")) +
    new_scale_color() +
    geom_tippoint(aes(size = iscow), color="black") +
    scale_size_manual(values=c(1, 2.5), name="cow isolate") +
    new_scale("size")+
    geom_tippoint(aes(color = iscow, size=iscow)) +
    scale_size_manual(values=c(0.5, 1.5), name="cow isolate") +
  
    scale_color_manual(values=c("grey60", "#238B45"), name="cow isolate") +
    theme(
      legend.position      = c(0, 0.95),   # (x, y) in npc units: 0 = left, 1 = top
      legend.justification = c(0, 0.95),   # anchor the legend’s own top-left corner
      legend.box.just      = "left",    # keep multiple guides left-aligned
      ## remove the white panels
      legend.background    = element_blank(),   # outer box
      legend.key           = element_blank(),   # keys behind symbols
      legend.box.background = element_blank()   # box around grouped legends
  
    ) +
    # label every year form 2020 to 2025
    scale_x_continuous(breaks=seq(2020, 2025, 1), labels=seq(2020, 2025, 1)) +
    ylim(1, 450) 
  
  p_tree <- p_tree +
    annotate("rect", xmin = 2020, xmax = 2020.5, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 1) +
    annotate("rect", xmin = 2021, xmax = 2021.5, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 1) +
    annotate("rect", xmin = 2022, xmax = 2022.5, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 1) +
    annotate("rect", xmin = 2023, xmax = 2023.5, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 1)+
    annotate("rect", xmin = 2024, xmax = 2024.5, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 1)+
    annotate("rect", xmin = 2025, xmax = 2025.5, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 1)
  
  # move the three newest layers (rectangles) to the bottom
  p_tree$layers <- append(p_tree$layers[ (length(p_tree$layers)-5) : length(p_tree$layers) ],
                          p_tree$layers[ 1 : (length(p_tree$layers)-6) ])
  

  if (ind){
    hpai_events = read.table("./combined/HLHxNx.constant.clades.tsv", header=TRUE, sep="\t")
  }else{
    hpai_events = read.table("./combined/HLHxNx.dependent.clades.tsv", header=TRUE, sep="\t")
  }
    
  # remove all events for which Lineage is not HPAI
  hpai_events$Time = mrsi - hpai_events$Height * 365.25
  # remove any event before 2020
  hpai_events = hpai_events[hpai_events$Time >= as.Date("2020-01-01"), ]
  
  lineage_type = c("HPAI", "LPAI")
  event_types = c("HPAI","HPAI+LPAI", "LPAI")
  distr = data.frame()
  for (i in seq(min(hpai_events$Sample), max(hpai_events$Sample))){
    for (l in lineage_type){
      # get the number of events for this sample and lineage
      n_events = hpai_events[hpai_events$Sample == i & hpai_events$Lineage == l, ]
      for (e in event_types){
        # check if lineages is in e, otherwise, next
        if (!grepl(l, e)){
          next
        }
        # get the number of events for this sample and lineage
        n_events_e = nrow(n_events[n_events$Event == e, ])
        # add to the distribution
        distr = rbind(distr, data.frame(
          Sample = i,
          Lineage = l,
          Event = e,
          n_events = n_events_e
        ))
      }
    }
  }
  
  p = ggplot(distr, aes(x=Event, y=n_events)) +
    geom_violin() +
    xlab("Event type") +
    ylab("Number of events") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    facet_grid(Lineage ~ ., scales = "free", switch = "y")
  plot(p)
  
  # for all events on HPAI lineages count how often each segment moved the other direction
  co_rea = data.frame()
  for (i in seq(min(hpai_events$Sample), max(hpai_events$Sample))){
    # get the number of events for this sample and lineage
    for (l in lineage_type){
    
      n_events = hpai_events[hpai_events$Sample == i & hpai_events$Lineage == l, ]
      for (e in event_types){
        # check if lineages is in e, otherwise, next
        if (!grepl(l, e)){
          next
        }
        # get the number of events for this sample and lineage
        n_events_e = n_events[n_events$Event == e, ]
        
        # loop over all segments not HA
        for (j in seq(2, length(segment_order))){
          # count how often i is in n_events_e$Segments
          count = sum(grepl(j-1, n_events_e$Segments))
          # add to the distribution
          co_rea = rbind(co_rea, data.frame(
            Sample = i,
            Lineage = l,
            Event = e,
            n_events = count,
            segment = segment_order[j]
          ))
        }
      }
    }
  }
  # calculate the quantiles for co_rea counts for each lineage, event and segment
  co_rea_quantiles = data.frame()
  for (s in unique(co_rea$segment)){
    for (e in unique(co_rea$Event)){
      for (l in unique(co_rea$Lineage)){
        # get the counts for this segment, event and lineage
        counts = co_rea[co_rea$segment == s & co_rea$Event == e & co_rea$Lineage == l, "n_events"]
        # get the quantiles
        q = 0.05
        lower = quantile(counts, q/2)
        upper = quantile(counts, 1-q/2)
        mean = mean(counts)
        co_rea_quantiles = rbind(co_rea_quantiles, data.frame(
          segment = s,
          Event = e,
          Lineage = l,
          mean = mean,
          lower = lower,
          upper = upper
        ))
      }
    }
  }
  
  lineage_colors <- c(
    "HPAI"    = "#E41A1C",  # red
    "LPAI"    = "#377EB8"  # blue
  )
  
  # convert Event HPAI to 'with other HPAI' and 'HPAI+LPAI' to 'coming from LPAI'
  co_rea_quantiles$evname = "HPAI"
  co_rea_quantiles$evname[co_rea_quantiles$Event == "HPAI+LPAI"] = "LPAI"
  
  # plot the posterior distribution of the number of events for each segment, facetting by event
  p_co = ggplot(co_rea_quantiles[co_rea_quantiles$Lineage=="HPAI", ], aes(x=segment, y=mean, color=evname)) +
    geom_point(position=position_dodge(-0.3)) +
    geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(-0.3), width=0.2) +
    ylab("")+
    scale_color_manual(values=lineage_colors, name="lineage origin of segment") +
    ylab("events") +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.box = "vertical"
    ) +
    guides(color = guide_legend(ncol = 1))  # Force legend to 1 column
    
  plot(p_co)
  
  # keeps track of teh moving average
  days_around = 90
  lineage_type = c("HPAI", "LPAI")
  event_types = c("HPAI","HPAI+LPAI", "LPAI")
  
  # moving_avg = data.frame()
  # for (l in lineage_type){
  #   # loop weekly from 2020 + 3 months to mrsi-3 months
  #   for (i in seq(as.Date("2020-04-01"), mrsi - days_around, by = "week")){
  #     # get the number of events for this sample and lineage
  #     n_events = hpai_events[hpai_events$Time >= i - days_around & 
  #                              hpai_events$Time <= i + days_around & 
  #                              hpai_events$Lineage == l, ]
  #     
  #     # for each sample, count how many events there are
  #     sum_events = c()
  #     for (s in unique(n_events$Sample)){
  #       # get the events for this sample
  #       sum_events = c(sum_events, nrow(n_events[n_events$Sample == s, ]))
  #     }
  #     
  #     
  #     
  #     # get the 95% HPD
  #     for (q in seq(0.05, 1.0, 0.05)){
  #       # get the quantiles
  #       moving_avg = rbind(moving_avg, data.frame(
  #         Time = i,
  #         Lineage = l,
  #         lower = quantile(sum_events, q/2),
  #         upper = quantile(sum_events, 1-q/2),
  #         q = q,
  #         type="All events involving HPAI offsprings"
  #       ))
  #     }
  #     
  #     for (e in event_types){
  #       # check if lineages is in e, otherwise, next
  #       if (!grepl(l, e)){
  #         next
  #       }
  #       # get the number of events for this sample and lineage
  #       n_events_e = n_events[n_events$Event == e, ]
  #       
  #       # for each sample, count how many events there are
  #       sum_events = c()
  #       for (s in unique(n_events$Sample)){
  #         # get the events for this sample
  #         n = nrow(n_events_e[n_events_e$Sample == s, ])
  #         # d = nrow(n_events[n_events$Sample == s, ])
  #         sum_events = c(sum_events, n)
  #       }
  #       
  #       
  #       # get the 95% HPD
  #       for (q in seq(0.05, 1.0, 0.05)){
  #         # get the quantiles
  #         moving_avg = rbind(moving_avg, data.frame(
  #           Time = i,
  #           Lineage = l,
  #           # lower = mean(sum_events),
  #           # upper = mean(sum_events),
  #           lower = quantile(sum_events, q/2),
  #           upper = quantile(sum_events, 1-q/2),
  #           q = q,
  #           type=e
  #         ))
  #       }
  #     }
  #   }
  # }
  
  # rename types to events between HPAI and LPAI and events only involving HPAI lineages
  # moving_avg$type[moving_avg$type == "HPAI+LPAI"] = "Events between HPAI and LPAI lineages"
  # moving_avg$type[moving_avg$type == "HPAI"] = "Events only involving HPAI lineages"
  # moving_avg$type[moving_avg$type == "LPAI"] = "Events only involving LPAI lineages"
  
  # moving_avg$Time = as.Date(moving_avg$Time)
  # 
  # dat = moving_avg[moving_avg$Lineage=="HPAI", ]
  # dat2 = moving_avg[moving_avg$Lineage=="HPAI" &
  #                     moving_avg$type == "Events between HPAI and LPAI lineages",]
  # 
  # 
  # # plot the moving average
  # p3 = ggplot(data = dat, aes(x=Time, y=upper, group=interaction(q, type))) +
  #   geom_ribbon(aes(ymin=lower, ymax=upper, fill=type), alpha=0.2) +
  #   scale_fill_manual(values=c("HPAI"="#E41A1C", "HPAI+LPAI"="#377EB8"), name="Type") +
  #   
  #   geom_line(data=dat2, aes(x=Time, y=upper*20), color="#377EB8", size=0.5) +
  #   scale_y_continuous(sec.axis = sec_axis(~./20, name="average proportion of\nevents with LPAI lineages")) +
  #   new_scale_color()+
  #   geom_line(data=smoothed_case_data, aes(x=date, y=positivity*20, color=type, group=type), method=NA, size=0.5) +
  #   scale_color_manual(values=c("HPAI"="#E41A1C", "LPAI"="#377EB8"), name="Type") +
  #   ylab("\n moving average of\nthe number of reassortment events") +
  #   coord_cartesian(ylim=c(0, 20)) +
  #   theme_minimal() 
  # plot(p3)
  
  # flip the axis of p_co
  
  
  
  # # plot the moving average
  # p = ggplot(moving_avg[moving_avg$Lineage=="HPAI", ], aes(x=Time, y=upper, group=q)) +
  #   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="#546E7A") +
  #   ylab("180 day moving averge of\nthe number of reassortment events\nwith HPAI offspring") +
  #   theme_minimal() + 
  #   new_scale_color() +
  #   ylim(0, 20) +
  #   facet_grid(Lineage~type)
  # plot(p)
  # 
  # ggsave(p, filename="../../Figures/h5n1_HPAI_moving_average.pdf", width=8, height=6)
  # # plot the moving average
  # p = ggplot(moving_avg[moving_avg$Lineage=="LPAI", ], aes(x=Time, y=upper, group=q)) +
  #   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="#546E7A") +
  #   ylab("180 day moving averge of\nthe number of reassortment events\nwith HPAI offspring") +
  #   theme_minimal() +
  #   facet_grid(Lineage~type)
  # plot(p)
  # 
  # ggsave(p, filename="../../Figures/h5n1_LPAI_moving_average.pdf", width=8, height=6)
  
  lineage_colors <- c(
    HPAI    = "#E41A1C",  # red
    LPAI    = "#377EB8",  # blue
    unknown  = "#4DAF4A"  # green
  )
  
  tree_plots = list()
  # plot all segments other than HA in the same plot
  for (s in seq(2, length(segment_order), 1)){
    # get the segment name
    segment = segment_order[s]
    # read in the tree
    tree = read.beast(paste0("./combined/HLHxNx.constant.", segment, ".tree"))
    # plot the tree, coloring by the trait HPAI+LPAI
    
    # convert mrsi to decimal year
    mrsi_dec <- as.numeric(mrsi - as.Date("1970-01-01")) / 365.25 + 1970
    
    
    # now plot
    p <- ggtree(tree, aes(color = type), mrsd = mrsi, size = 0.2) +
      theme_tree2() +
      coord_cartesian(xlim=c(2020, mrsi_dec+0.05))+
      # add cow tips as filled circles
      scale_color_manual(values=lineage_colors, name="Pathogenicity\nbased on HA", na.translate  = FALSE)+
      theme(
        legend.position      = c(0, 1),   # (x, y) in npc units: 0 = left, 1 = top
        legend.justification = c(0, 1),   # anchor the legend’s own top-left corner
        legend.box.just      = "left",    # keep multiple guides left-aligned
        ## remove the white panels
        legend.background    = element_blank(),   # outer box
        legend.key           = element_blank(),   # keys behind symbols
        legend.box.background = element_blank()   # box around grouped legends
        
      ) +
      scale_x_continuous(breaks = seq(2020, 2025, 2), labels = seq(2020, 2025, 2))+  # Every 2 years
      # new_scale_color() +
      geom_tippoint(color="black", size=1) +
      # scale_size_manual(values=c(1, 2.5), name="cow isolate") +
      # new_scale("size")+
      geom_tippoint(aes(color = type), size=0.5)+
      ylim(1, 450) 

    
    p <- p +
      annotate("rect", xmin = 2020, xmax = 2020.5, ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 1) +
      annotate("rect", xmin = 2021, xmax = 2021.5, ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 1) +
      annotate("rect", xmin = 2022, xmax = 2022.5, ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 1) +
      annotate("rect", xmin = 2023, xmax = 2023.5, ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 1)+
      annotate("rect", xmin = 2024, xmax = 2024.5, ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 1)+
      annotate("rect", xmin = 2025, xmax = 2025.5, ymin = -Inf, ymax = Inf,
               fill = "grey95", alpha = 1)
    
    # move the three newest layers (rectangles) to the bottom
    p$layers <- append(p$layers[ (length(p$layers)-5) : length(p$layers) ],
                       p$layers[ 1 : (length(p$layers)-6) ])
    
    tree_plots[[s-1]] <- p
  }
  
  p_righ <- plot_grid(tree_plots[[5-1]]+theme(legend.position="none"), tree_plots[[7-1]]+theme(legend.position="none"), labels=c('C','D'), 
                      ncol=2, align="h", axis="l")
  p_top <- plot_grid(p_co, p_righ, labels=c('B',''), ncol=1, align="h", axis="l", rel_heights=c(0.5, 1))
  pcomp <- plot_grid(p_tree, p_top, labels=c('A',''), ncol=2, align="h", axis="l", rel_widths=c(1, 0.5))
  plot(pcomp)
  
  if (ind){
    ggsave(pcomp, filename="../../Figures/Figure2.pdf", width=12, height=8)
  }else{
    ggsave(pcomp, filename="../../Figures/h5n1_reassortment_dependent.pdf", width=12, height=8)
  }
  
  
  
  # 1) extract the legend from one tree plot
  legend_grob <- get_legend(
    tree_plots[[1]] +
      theme(
        legend.position      = "right",
        legend.justification = "center",
        legend.key.size      = unit(2, "lines"),          # make the keys taller
        legend.title         = element_text(size = 16),    # bigger title
        legend.text          = element_text(size = 14),    # bigger labels
        legend.spacing.y     = unit(0.5, "lines")         # more spacing between entries
        
      )
  )
  
  # 2) turn it into a “plot” object
  legend_panel <- ggdraw(legend_grob)
  
  # remove legend from all tree_plots
  tree_plots <- lapply(tree_plots, function(p) {
    p + theme(legend.position = "none")  # remove legend from each tree plot
  })
  
  # 3) combine your 7 trees + the legend as the 8th panel
  all_panels <- c(tree_plots, list(legend_panel))
  
  # 4) arrange in a 2×4 grid
  combined_trees <- plot_grid(
    plotlist = all_panels,
    ncol     = 4,
    labels   = segment_order[-1],  # will label first 7 panels NA, MP, …, PA
    label_size = 12,
    align    = "hv"
  )
  
  # 5) render and save
  print(combined_trees)
  if (ind){
    ggsave(combined_trees,filename = "../../Figures/h5n1_all_segment_trees_constant.pdf",
           width = 12, height = 8)
  }else{
    ggsave(combined_trees,filename = "../../Figures/h5n1_all_segment_trees_dependent.pdf",
           width = 12, height = 8)
  }
}



