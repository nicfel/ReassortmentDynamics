
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
    TRUE ~ location))




tr_p <- read.beast("MCC/HPAI_HA_northamerica_targeted_dta.mcc.tree")

tip_data_p <- tibble(label = tr_p@phylo$tip.label) %>%
  mutate(strain = str_split_fixed(label, "\\|", 2)[,1])

treedata <- left_join(tip_data_p,metadata, by = "strain" )
genoflu_over10 <- names(which(table(treedata$genoFLU) > 10))



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
  filter(genoFLU %in% genoflu_over10) %>%
  filter(!str_detect(genoFLU, "Minor")) %>%
  filter(!str_detect(genoFLU, "Unseen")) %>%
  filter(!str_detect(genoFLU, "divergent")) %>%
  mutate(tree_file = str_remove_all(tree_file,"_northamerica_targeted_dta\\.mcc\\.tree|HPAI_")) %>%
  mutate(
    location = ifelse(
      location == "mississippi_flyway+central_flyway+",
      "mississippi_flyway",
      location),
    tree_file = factor(tree_file, levels = segment_order))

flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F")

p <- ggplot()
locations <- unique(plot_df$location)

for (i in seq_along(locations)) {
  loc <- locations[i]
  df_loc <- plot_df %>% filter(location == loc)
  
  if (i > 1) p <- p + new_scale("alpha")
  
  p <- p +
    geom_errorbarh(
      data = df_loc,
      aes(
        y = tree_file,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high,
        alpha = location_prob),
      color = flyway_colors[loc],
      height = 0.25,
      linewidth = 1) +
    geom_point(
      data = df_loc,
      aes(
        x = tmrca_decimal_date,
        y = tree_file,
        alpha = location_prob),
      color = flyway_colors[loc],
      size = 3) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
  
  df_grey <- df_loc %>% filter(location_prob < 0.5)
  
  p <- p +
    geom_errorbarh(data = df_grey, aes(
        y = tree_file,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high),
      color = "grey70",
      height = 0.25,
      linewidth = 1) +
    geom_point(data = df_grey, aes(x = tmrca_decimal_date,
        y = tree_file),
      color = "grey70",
      size = 3)
}

pf <- p +
  facet_wrap(~ genoFLU, ncol = 4) +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = "Date",
    y = "Segment") +
  theme_minimal(base_size = 14)


pf

ggsave("flyway_combinedtmcra_sub10_pp50-facet.pdf", plot = pf, height = 12, width = 14, units = "in")

##################
segment_order <- c("NS", "MP", "NA", "NP", "HA", "PA", "PB1", "PB2")

plot_df <- all_results %>%
  filter(genoFLU %in% genoflu_over10) %>%
  filter(!str_detect(genoFLU, "Minor")) %>%
  filter(!str_detect(genoFLU, "Unseen")) %>%
  filter(!str_detect(genoFLU, "divergent")) %>%
  mutate(tree_file = str_remove_all(tree_file, "_northamerica_targeted_dta\\.mcc\\.tree|HPAI_")) %>%
  mutate(location = ifelse(location == "central_flyway+mississippi_flyway","mississippi_flyway",location),
    tree_file = factor(tree_file, levels = segment_order),
    genoFLU = factor(genoFLU, levels = sort(unique(genoFLU))),
    geno_segment = paste(genoFLU, tree_file, sep = " | "))

flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F")

p <- ggplot()
locations <- unique(plot_df$location)

for (i in seq_along(locations)) {
  loc <- locations[i]
  df_loc <- plot_df %>% filter(location == loc)
  
  if (i > 1) p <- p + new_scale("alpha")
  
  p <- p +
    geom_errorbarh(data = df_loc,aes(
        y = geno_segment,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high,
        alpha = location_prob),
      color = flyway_colors[loc],
      height = 0.2,
      linewidth = 1.2) +
    geom_point(
      data = df_loc,
      aes(
        x = tmrca_decimal_date,
        y = geno_segment,
        alpha = location_prob),
      color = flyway_colors[loc],
      size = 2) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
  
  df_grey <- df_loc %>% filter(location_prob < 0.5)
  
  p <- p +
    geom_errorbarh(data = df_grey,
      aes(
        y = geno_segment,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high),
      color = "grey70",
      height = 0.2,
      linewidth = 1.2) +
    geom_point(
      data = df_grey,
      aes(
        x = tmrca_decimal_date,
        y = geno_segment),
      color = "grey70",
      size = 2)
}

p1 <- p +
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

p1

ggsave("flyway_combinedtmcra_sub10_pp50.pdf", plot = p1, height = 8.5, width = 12, units = "in")

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
