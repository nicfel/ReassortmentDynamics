
library(treeio)
library(tidytree)
library(ape)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(patchwork)
library(viridis)


setwd("~/Dropbox/targeted-beast-DTA-HPAI-reassortment/h5-data-updates-main-reassortmentprev/h5nx/")
most_recent_date <- 2025.33   # 2025-05-01 mrsd
tree_folder <- "MCC/"
metadata_file <- "metadata-with-clade.tsv"
# metadta file comes from moncla lab github repo for nextstrain

# Read metadata
metadata <- read_tsv(metadata_file, col_types = cols()) %>%
  rename(strain = strain, genoFLU = genoflu)

# Function to calculate TMRCA per genoFLU group
process_tree <- function(tree_file, metadata, most_recent_date) {
  cat("Processing:", tree_file, "\n")
  tr <- read.beast(tree_file)
  
  # Extract tip labels and clean strain names
  tip_data <- tibble(label = tr@phylo$tip.label) %>%
    mutate(strain = str_split_fixed(label, "\\|", 2)[,1])
  
  # Join with metadata
  tip_data <- tip_data %>%
    left_join(metadata, by = "strain")
  
  # Filter out missing genoFLU
  present <- tip_data %>% filter(!is.na(genoFLU))
  
  # For each genoFLU group, find MRCA node number
  res <- present %>%
    group_by(genoFLU) %>%
    summarise(
      mrca_node = MRCA(tr@phylo, label),
      .groups = "drop") %>%
    # make sure MRCA is numeric (node index in phylo)
    mutate(mrca_node = as.integer(mrca_node)) %>%
    # join with tr@data (node annotation table)
    left_join(tr@data %>%
                mutate(node = as.integer(node)),
              by = c("mrca_node" = "node")) %>%
    mutate(
      tmrca_height_years_bp = as.numeric(height),
      tmrca_decimal_date = most_recent_date - tmrca_height_years_bp,
      height_range = as.character(height_range),
      height_range = str_remove_all(height_range, "[c()]"),
      height_low = as.numeric(str_split_fixed(height_range, ",", 2)[,1]),
      height_high = as.numeric(str_split_fixed(height_range, ",", 2)[,2]),
      
      # decimal dates for bounds
      tmrca_decimal_date_low = most_recent_date - height_low,
      tmrca_decimal_date_high = most_recent_date - height_high,
      location = location,
      location_prob = as.numeric(location.prob)) %>%
    select(genoFLU, tmrca_height_years_bp, tmrca_decimal_date, height_range, height_low, height_high, 
           tmrca_decimal_date_low, tmrca_decimal_date_high, location, location_prob)
  
  res$tree_file <- basename(tree_file)
  return(res)
}

# Process all trees in folder
tree_files <- list.files(tree_folder, pattern = "*.tree$", full.names = TRUE)

all_results <- map_dfr(tree_files, ~ process_tree(.x, metadata, most_recent_date))


write_csv(all_results, "tmrca_results.csv")

all_results
all_results <- all_results %>%
  mutate(location = case_when(
    location == "atl" ~ "atlantic_flyway",
    location == "cen" ~ "central_flyway",
    location == "mis" ~ "mississippi_flyway",
    location == "pac" ~ "pacific_flyway",
    location == "mis+cen" ~ "mississippi_flyway",
    TRUE ~ location
  ))


# Function to calculate MRCA-to-first-taxon branch length per genoFLU group
process_tree_first_taxon <- function(tree_file, metadata, most_recent_date) {
  cat("Processing (first taxon):", tree_file, "\n")
  tr <- read.beast(tree_file)
  
  # Extract tip labels and clean strain names
  tip_data <- tibble(label = tr@phylo$tip.label) %>%
    mutate(strain = str_split_fixed(label, "\\|", 2)[,1])
  
  # Join with metadata
  tip_data <- tip_data %>%
    left_join(metadata, by = "strain")
  
  # Filter out missing genoFLU
  present <- tip_data %>% filter(!is.na(genoFLU))
  
  # Precompute distance matrix
  dmat <- dist.nodes(tr@phylo)
  
  # For each genoFLU group, compute MRCA → closest tip distance
  res <- present %>%
    group_by(genoFLU) %>%
    summarise(
      mrca_node = MRCA(tr@phylo, label),
      .groups = "drop") %>%
    mutate(mrca_node = as.integer(mrca_node)) %>%
    rowwise() %>%
    mutate(
      group_tip_indices = list(match(present$label[present$genoFLU == genoFLU], tr@phylo$tip.label)),
      first_taxon_rel_years = min(dmat[mrca_node, unlist(group_tip_indices)])) %>%
    ungroup() %>%
    select(genoFLU, first_taxon_rel_years)
  
  res$tree_file <- basename(tree_file)
  return(res)
}



all_first_taxon_results <- map_dfr(tree_files, ~ process_tree_first_taxon(.x, metadata, most_recent_date))

write_csv(all_first_taxon_results, "first_taxon_results.csv")
all_first_taxon_results




###########
#pltos for time from mcra to first detected sample. 


plot_df <- all_first_taxon_results %>%
  mutate(first_taxon_rel_days = first_taxon_rel_years * 365.25) %>%
  filter(
    genoFLU != "Unseen constellation",
    genoFLU != "Not assigned (too divergent)",
    !str_detect(genoFLU, "Minor"))

plot_df <- plot_df %>%
  mutate(tree_file = str_remove_all(tree_file, "_northamerica_targeted_dta\\.mcc\\.tree|HPAI_"))


# all together 
ggplot(plot_df, aes(x = first_taxon_rel_days)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30, color = "black", alpha = 0.6) +
  geom_density(color = "darkgreen", size = 1.2, ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "MRCA → first taxon distance (days)",
    y = "Density",
    title = "Distribution of MRCA-to-first-taxon waiting times")


##### facted
ggplot(plot_df, aes(x = first_taxon_rel_days)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30, color = "black", alpha = 0.6) +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600), minor_breaks = waiver()) +
  geom_density(color = "darkgreen", size = 1.2) +
  theme_minimal(base_size = 18) +
  labs(
    x = "Days",
    y = "Density",
    title = "Distribution of time to first detection") +
  facet_wrap(~ tree_file, scales = "fixed")




ggplot(plot_df, aes(x = first_taxon_rel_days, color = tree_file, fill = tree_file)) +
  geom_density(alpha = 0.2, size = 1.2) +
  theme_minimal(base_size = 18) +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600), minor_breaks = waiver()) +
  labs(
    x = "Days",
    y = "Density",
    title = "Distribution of time to first detection",
    color = "Segment",
    fill  = "Segment") +
  scale_color_viridis(discrete = TRUE, option = "C", direction = -1) 





##############
#Visualization of the results for single genotype
target_geno <- "B3.13"
plot_df <- all_results %>%
  filter(genoFLU == target_geno)

plot_df <- plot_df %>%
  mutate(tree_file = str_remove_all(tree_file, "_northamerica_targeted_dta\\.mcc\\.tree|HPAI_"))


flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F"
)

locations <- unique(plot_df$location)

p <- ggplot()

for (i in seq_along(locations)) {
  loc <- locations[i]
  df_loc <- plot_df %>% filter(location == loc)
  
  if (i > 1) p <- p + new_scale("alpha")
  
  p <- p +
    geom_point(
      data = df_loc,
      aes(x = tmrca_decimal_date,y = tree_file,alpha = location_prob),
      color = flyway_colors[loc],
      size = 4) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
}

p +
  labs(title = paste("TMRCA estimates for", target_geno),
  x = "Date",
  y = "Tree file") +
  theme_minimal(base_size = 14)



####################################
segment_order <- c("NS", "MP", "NA", "NP", "HA", "PA", "PB1", "PB2")


plot_df <- all_results %>%
  filter(!str_detect(genoFLU, "Minor")) %>%
  filter(!str_detect(genoFLU, "Unseen")) %>%
  filter(!str_detect(genoFLU, "divergent")) %>%
  mutate(tree_file = str_remove_all(tree_file, "_northamerica_targeted_dta\\.mcc\\.tree|HPAI_")) %>%
  mutate(location = ifelse(location == "mississippi_flyway+central_flyway+", "mississippi_flyway", location))

plot_df <- plot_df %>%
  mutate(
tree_file = factor(tree_file, levels = segment_order))
                   
                   

flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F"
)

# build layers across all locations
p <- ggplot()
locations <- unique(plot_df$location)

for (i in seq_along(locations)) {
  loc <- locations[i]
  df_loc <- plot_df %>% filter(location == loc) 
  
  if (i > 1) p <- p + new_scale("alpha")
  
  
  p <- p +
    geom_point(
      data = df_loc,
      aes(
        x = tmrca_decimal_date,
        y = tree_file, alpha = location_prob), color = flyway_colors[loc], size = 5) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
}

p
  
p +
  facet_wrap(~ genoFLU, ncol = 4 ) +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = "Date",
    y = "Segement") + 
  theme_minimal(base_size = 14) 

########################
plot_df <- all_results %>%
  filter(!str_detect(genoFLU, "Minor")) %>%
  filter(!str_detect(genoFLU, "Unseen")) %>%
  filter(!str_detect(genoFLU, "divergent")) %>%
  mutate(tree_file = str_remove_all(tree_file, "_northamerica_targeted_dta\\.mcc\\.tree|HPAI_")) %>%
  mutate(
    genoFLU = factor(genoFLU, levels = sort(unique(genoFLU))),
    geno_segment = paste(genoFLU, tree_file, sep = " | ")) %>%
  mutate(location = ifelse(location == "central_flyway+mississippi_flyway", "mississippi_flyway", location))


x_range <- max(plot_df$tmrca_decimal_date, na.rm = TRUE) - min(plot_df$tmrca_decimal_date, na.rm = TRUE)
x_offset <- ifelse(is.finite(x_range) & x_range > 0, x_range * 0.02, 0.01)

# Build HA label table: only genotypes that have a HA row will get a label
ha_labels <- plot_df %>%
  filter(str_detect(tree_file, regex("HA$", ignore_case = TRUE))) %>%
  group_by(genoFLU) %>%
  ungroup() %>%
  mutate(
    x_label = tmrca_decimal_date + x_offset,
    y_label = paste(genoFLU, tree_file, sep = " | ")
  )

# plotting

flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F"
)

p <- ggplot()
locations <- unique(plot_df$location)

for (i in seq_along(locations)) {
  loc <- locations[i]
  df_loc <- plot_df %>% filter(location == loc)
  
  if (i > 1) p <- p + ggnewscale::new_scale("alpha")
  
  p <- p +
    geom_errorbarh(
      data = df_loc,
      aes(
        y = geno_segment,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high,
        alpha = location_prob),
      color = flyway_colors[loc],
      height = 0.2,
      linewidth = 1.2,
      show.legend = FALSE) +
    geom_point(
      data = df_loc,
      aes(
        x = tmrca_decimal_date,
        y = geno_segment,
        alpha = location_prob),
      color = flyway_colors[loc],
      size = 4,
      show.legend = FALSE) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
}

# add genotype labels next to each genotype's HA point
p <- p +
  geom_text(
    data = ha_labels,
    aes(x = x_label, y = y_label, label = genoFLU),
    inherit.aes = FALSE,
    hjust = 0,
    fontface = "bold",
    size = 5) +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = "Date",
    y = "Genotype | Segment") +
  theme_minimal(base_size = 20) +
  theme(axis.text.y = element_text(size = 8))

p


##################
segment_order <- c("NS", "MP", "NA", "NP", "HA", "PA", "PB1", "PB2")


plot_df <- all_results %>%
  filter(!str_detect(genoFLU, "Minor")) %>%
  filter(!str_detect(genoFLU, "Unseen")) %>%
  filter(!str_detect(genoFLU, "divergent")) %>%
  mutate(tree_file = str_remove_all(tree_file, "_northamerica_targeted_dta\\.mcc\\.tree|HPAI_")) %>%
  mutate(location = ifelse(location == "central_flyway+mississippi_flyway", "mississippi_flyway", location))


plot_df <- plot_df %>%
  mutate(
    tree_file = factor(tree_file, levels = segment_order))


flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F"
)

plot_df <- plot_df %>%
  mutate(
    tree_file = factor(tree_file, levels = segment_order),
    genoFLU = factor(genoFLU, levels = sort(unique(genoFLU))),
    geno_segment = paste(genoFLU, tree_file, sep = " | "))


p <- ggplot()
locations <- unique(plot_df$location)

for (i in seq_along(locations)) {
  loc <- locations[i]
  df_loc <- plot_df %>% filter(location == loc)
  
  if (i > 1) p <- p + new_scale("alpha")
  
  p <- p +
    # horizontal error bars
    geom_errorbarh(
      data = df_loc,
      aes(
        y = geno_segment,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high,
        alpha = location_prob),
      color = flyway_colors[loc],
      height = 0.2,
      linewidth = 1.2) +
    # points
    geom_point(
      data = df_loc,
      aes(
        x = tmrca_decimal_date,
        y = geno_segment,
        alpha = location_prob),
      color = flyway_colors[loc],
      size = 2) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
}


p +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = "Date",
    y = NULL) +
  facet_grid(
    genoFLU ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"),
    strip.background = element_blank(),
    panel.spacing.y = unit(0.025, "lines")) +
  scale_x_continuous(limits = c(2019, 2025))


p_top <- p +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = NULL,
    y = NULL) +
  facet_grid(
    genoFLU ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"),
    strip.background = element_blank(),
    panel.spacing.y = unit(0.0000005, "lines"),
    legend.position = c(0.05, 0.05),          
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key.size = unit(0.4, "cm"),  
    legend.text = element_text(size = 8),     
    legend.title = element_text(size = 9)) +
  scale_x_continuous(limits = c(2020, 2025))

p_top
# Bottom plot (weekly counts)

# flywy colors defined
flyway_colors <- c(
  "Pacific Flyway" = "#E56C2F",
  "Central Flyway" = "#CEB541",
  "Mississippi Flyway" = "#69B091",
  "Atlantic Flyway" = "#4274CE"
)

# for dataframe for botoom panel with detections see detections_for_figure.R script

p_bottom <- ggplot(flyway_weekly_counts_long, aes(x = week, y = count, color = flyway)) +
  geom_line() +
  scale_color_manual(values = flyway_colors, name = "Flyway") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
    limits = as.Date(c("2020-01-01", "2025-01-01"))) +
  labs(
    x = "Week",
    y = "Number of Detections") +
  facet_wrap(~ path_type, ncol = 1, scales = "free_y") +
  ylim(0,200) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.05, 0.05),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12))

p_combined <- p_top / p_bottom +
  plot_layout(heights = c(4, 1))

p_combined

#######################################


##############################
# edge colorings based on the state
# Based on reweards  visualization
# use the MCC tree with edges annotated to the branches

beast_tree <- read.beast("MCC/ordgooseduck_MKJ_HA.mcc.tree")
data<-as_tibble(beast_tree)
## filter data
df<- data %>% select(ord,height)

# input date of earliest seq
df <- df %>% mutate(year = 2025.33 - as.numeric(height)) %>%
  mutate(ord = sub("\\+.*", "", ord))
## get the propotrion of rewards

interval_size <- 1/52

# Round or floor to weekly decimal intervals
df_long_weekly <- df %>%
  mutate(week_decimal = floor(year / interval_size) * interval_size) %>%
  group_by(week_decimal, ord) %>%
  summarise(total_freq = sum(as.numeric(height)), .groups = "drop")

all_weeks <- seq(min(df_long_weekly$week_decimal),
                 max(df_long_weekly$week_decimal),
                 by = interval_size)

all_ords <- sort(unique(df_long_weekly$ord))

df_long_weekly_complete <- df_long_weekly %>%
  complete(week_decimal = all_weeks, ord = all_ords, fill = list(total_freq = 0))

df_long_weekly_percent <- df_long_weekly_complete %>%
  group_by(week_decimal) %>%
  mutate(percent = total_freq / sum(total_freq) * 100) %>%
  ungroup()

p_week_decimal <- ggplot(df_long_weekly_percent, aes(x = week_decimal, y = percent, fill = ord)) +
  geom_area(size = 0.2, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0,0)) +
  xlim(2022,2025.33) +
  labs(
    title = "% Edges by Host Order (Weekly)",
    x = "Year",
    y = "Proportion",
    fill = "Order") +
  theme_minimal(base_size = 18) +
  scale_fill_ptol()

p_week_decimal

# Save the plot
ggsave("stacked_week_order.pdf", plot = p_week_decimal, height = 8.5, width = 11, units = "in")



####################

all_results <- read.csv("tmrca_results_flyway.csv")
mixed <- all_results %>%
  filter(!str_detect(genoFLU, "Minor")) %>%
  filter(!str_detect(genoFLU, "Unseen")) %>%
  filter(!str_detect(genoFLU, "divergent")) %>%
  group_by(genoFLU) %>%
  filter(n_distinct(location) > 1) %>%
  ungroup()

combo_df <- mixed %>%
  group_by(genoFLU) %>%
  summarise(
    loc_combo = location 
    %>% unique() 
    %>% sort() 
    %>% paste(collapse = "+"),.groups = "drop")

loc_combo_proportions <- combo_df %>%
  count(loc_combo, name = "count") %>%
  mutate(proportion = count / sum(count))

loc_combo_proportions

write.csv(loc_combo_proportions,"loc_combo_genoflu_proportions.csv", row.names = FALSE)

