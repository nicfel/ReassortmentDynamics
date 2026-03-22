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
library(scico)
library(ggthemes)
library(scales)
library(zoo)
# script to process MCC trees for genotype flyways analysis 

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


### facted
ggplot(plot_df, aes(x = first_taxon_rel_days)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30, color = "black", alpha = 0.6) +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600), minor_breaks = waiver()) +
  geom_density(color = "darkgreen", size = 1.2) +
  theme_minimal(base_size = 18) +
  labs(x = "Days",y = "Density",title = "Distribution of time to first detection") +
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
    tree_file = factor(tree_file, levels = segment_order),
    genoflu = factor(genoFLU, levels = sort(unique(genoFLU))),
    geno_segment = paste(genoFLU, tree_file, sep = " | "))

write.csv(plot_df,"plot_df_flyway.csv", row.names = FALSE)

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
    geom_errorbarh(data = df_loc,aes(y = tree_file,xmin = tmrca_decimal_date_low,xmax = tmrca_decimal_date_high,alpha = location_prob),color = flyway_colors[loc],height = 0.25,linewidth = 1) +
    geom_point(data = df_loc,aes(x = tmrca_decimal_date,y = tree_file,alpha = location_prob),color = flyway_colors[loc],size = 3) +
    scale_alpha(range = c(0.3, 1), name = paste0(loc, " prob"))
  
  df_grey <- df_loc %>% filter(location_prob < 0.5)
  
  p <- p +
    geom_errorbarh(data = df_grey, aes(y = tree_file,xmin = tmrca_decimal_date_low,xmax = tmrca_decimal_date_high),
      color = "grey70",
      height = 0.25,
      linewidth = 1) +
    geom_point(data = df_grey, aes(x = tmrca_decimal_date,y = tree_file), color = "grey70", size = 3)
}

pf <- p +
  facet_wrap(~ genoFLU, ncol = 4) +
  labs(title = "TMRCA estimates across all genotypes (excluding Minor)",x = "Date",y = "Segment") +
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
  labs(title = "TMRCA estimates across all genotypes (excluding Minor)",x = "Date",y = NULL) +
  facet_grid(genoFLU ~ .,scales = "free_y",space = "free_y",switch = "y") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),strip.placement = "outside",strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"),strip.background = element_blank(),panel.spacing.y = unit(0.025, "lines")) +
  scale_x_continuous(limits = c(2019, 2025))

p1

ggsave("flyway_combinedtmcra_sub10_pp50.pdf", plot = p1, height = 8.5, width = 12, units = "in")



#####

beast_tree <- read.beast("MCC/HPAI_HA_northamerica_targeted_dta.mcc.tree")
data<-as_tibble(beast_tree)
data <- data %>% 
  mutate(location = sub("\\+.*", "", location))
## filter data
df<- data %>% select(location,height)

# input date of earliest seq
df <- df %>% mutate(year = 2025.33 - as.numeric(height)) %>%
  mutate(location = sub("\\+.*", "", location))
## get the propotrion of rewards

interval_size <- 1/52

df_long_weekly <- df %>%
  mutate(
    week_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  group_by(week_decimal, location) %>%
  summarise(total_freq = n(), .groups = "drop")

all_weeks <- sort(unique(df_long_weekly$week_decimal))
all_locations  <- sort(unique(df_long_weekly$location))

df_long_weekly_complete <- df_long_weekly %>%
  complete(week_decimal = all_weeks, location = all_locations, fill = list(total_freq = 0))

df_long_weekly_percent <- df_long_weekly_complete %>%
  group_by(week_decimal) %>%
  mutate(
    prop = if (sum(total_freq) == 0) 0 else total_freq / sum(total_freq)) %>%
  ungroup()



flyway_colors <- c(
  atl    = "#4274CE",
  cen     = "#CEB541",
  mis = "#69B091",
  pac     = "#E56C2F")

p_week_decimal <- ggplot(df_long_weekly_percent, aes(x = week_decimal, y = prop, fill = location)) +
  geom_area(size = 0.2, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0,0)) +
  xlim(2021.5, 2025.33) +
  labs(title = "% Edges by Flyway (Weekly)",x = "Year",y = "Proportion",fill = "flyway") +
  theme_minimal(base_size = 18) +
  scale_fill_manual(values = flyway_colors)

p_week_decimal

ggsave("stacked_week_flyway.pdf", plot = p_week_decimal, height = 8.5, width = 12, units = "in")


#################################
interval_size <- 1 / 12
x_limits <- c(2021.5, 2025)
area_scale <- 1


geno_monthly <- plot_df %>%
  mutate(month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, name = "geno_n")

# panel A - donors of segements

location_prop <- plot_df %>%
  mutate(
    location = sub("\\+.*", "", location),
    month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, location, name = "n_location") %>%
  complete(
    month_decimal = sort(unique(month_decimal)),
    location = sort(unique(location)),
    fill = list(n_location = 0)) %>%
  group_by(month_decimal) %>%
  mutate(prop = if (sum(n_location) == 0) 0 else n_location / sum(n_location)) %>%
  ungroup()

location_scaled <- location_prop %>%
  left_join(geno_monthly, by = "month_decimal") %>%
  mutate(geno_n = replace_na(geno_n, 0),scaled_height = prop * geno_n * area_scale)




flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F")

p_panel_A <- ggplot() +
  geom_area(data = location_scaled,aes(x = month_decimal, y = scaled_height, fill = location),alpha = 0.9,linewidth = 0.2) +
  geom_line(data = geno_monthly,aes(x = month_decimal, y = geno_n),linewidth = 1.3,color = "black") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = flyway_colors) +
  labs(title = "A) Donor segment proportion",y = "Genotype segements",x = NULL,fill = "flyway") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20))

p_panel_A


# panel b - edges of the phylogeny 
edge_prop <- as_tibble(beast_tree) %>%
  mutate(
    location = sub("\\+.*", "", location),
    year = 2025.33 - as.numeric(height),
    month_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, location, name = "n_edges")



all_months <- seq(min(edge_prop$month_decimal), max(edge_prop$month_decimal), by = interval_size)
all_locations <- sort(unique(edge_prop$location))

edge_prop_complete <- edge_prop %>%
  complete(
    month_decimal = all_months,
    location = all_locations,
    fill = list(n_edges = 0)) %>%
  group_by(month_decimal) %>%
  mutate(prop = if (sum(n_edges) == 0) 0 else n_edges / sum(n_edges)) %>%
  ungroup()

edge_scaled <- edge_prop_complete %>%
  left_join(geno_monthly, by = "month_decimal") %>%
  arrange(month_decimal) %>%
  group_by(location) %>%
  mutate(geno_n = na.approx(geno_n, x = month_decimal, na.rm = FALSE, rule = 2),
         scaled_height = prop * geno_n * area_scale) %>%
  ungroup()



flyway_colors <- c(
  atl    = "#4274CE",
  cen     = "#CEB541",
  mis = "#69B091",
  pac     = "#E56C2F")


p_panel_B <- ggplot() +
  geom_area(data = edge_scaled,aes(x = month_decimal, y = scaled_height, fill = location),alpha = 0.9,linewidth = 0.2) +
  geom_line(data = geno_monthly,aes(x = month_decimal, y = geno_n),linewidth = 1.3,color = "black") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = flyway_colors) +
  labs(title = "B) Edges proportion",y = "Genotype segments",x = "Year",fill = "flyway") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20))


p_panel_B


final_figure <- p_panel_A / p_panel_B +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

final_figure

ggsave( "geno_edges_vs_donors_AB_flyways.pdf", plot = final_figure, height = 12, width = 12, units = "in")




###############################
# PANEL A — donor segments

flyway_colors <- c(
  atlantic_flyway    = "#4274CE",
  central_flyway     = "#CEB541",
  mississippi_flyway = "#69B091",
  pacific_flyway     = "#E56C2F")

p_panel_A <- ggplot() +
  geom_col(data = location_scaled,aes(x = month_decimal, y = scaled_height, fill = location),width = interval_size * 0.9)  +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = flyway_colors) +
  labs(title = "A) Donor segment proportion",y = "Genotype segments",x = NULL,fill = "flyway") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20))

##########
# PANEL B — phylogeny edges
flyway_colors <- c(
  atl    = "#4274CE",
  cen     = "#CEB541",
  mis = "#69B091",
  pac     = "#E56C2F")



p_panel_B <- ggplot() +
  geom_col(data = edge_scaled,aes(x = month_decimal, y = scaled_height, fill = location),width = interval_size * 0.9) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = flyway_colors) +
  labs(title = "B) Edges proportion",y = "Genotype segments",x = "Year",fill = "flyway") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20))

#################################
# combine panels

final_figure <- p_panel_A / p_panel_B +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

final_figure

ggsave("geno_edges_vs_donors_AB__bar_flyways.pdf",plot = final_figure,height = 12,width = 12,units = "in")


###########
# correlation calculation 

interval_size <- 1 / 12
x_limits <- c(2021.5, 2025.33)
cor_method <- "spearman"

# donor segments calc
location_prop_obs <- plot_df %>%
  mutate(location = sub("\\+.*", "", location),
    location = recode(location,
      atlantic_flyway    = "atl",
      mississippi_flyway = "mis",
      central_flyway     = "cen",
      pacific_flyway     = "pac"),
    month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, location, name = "n_obs") %>%
  group_by(month_decimal) %>%
  mutate(prop_obs = n_obs / sum(n_obs)) %>%
  ungroup() %>%
  select(month_decimal, location, prop_obs)



# edge proportions
location_prop_edge <- as_tibble(beast_tree) %>%
  mutate(
    location = sub("\\+.*", "", location),
    year = 2025.33 - as.numeric(height),
    month_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, location, name = "n_edge") %>%
  group_by(month_decimal) %>%
  mutate(prop_edge = n_edge / sum(n_edge)) %>%
  ungroup() %>%
  select(month_decimal, location, prop_edge)


# format data
location_long <- location_prop_obs %>%
  full_join(location_prop_edge,
            by = c("month_decimal", "location")) %>%
  replace_na(list(prop_obs = 0, prop_edge = 0)) %>%
  pivot_longer(cols = c(prop_obs, prop_edge),names_to = "source",values_to = "prop") %>%
  mutate(source = recode(source, prop_obs = "Donor", prop_edge = "Edges"))


# calc spearman corrl
location_cor <- location_prop_obs %>%
  inner_join(location_prop_edge,
             by = c("month_decimal", "location")) %>%
  group_by(location) %>%
  summarise(
    R = if (n() >= 1)
      cor(prop_obs, prop_edge, method = cor_method), .groups = "drop") %>%
  mutate(label = paste0("R = ", round(R, 2)))


p_faceted <- ggplot(
  location_long,
  aes(x = month_decimal, y = prop, color = source)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  facet_wrap(~ location, scales = "fixed") +
  scale_x_continuous(limits = x_limits,expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),expand = c(0, 0)) +
  scale_color_manual(values = c("Donor" = "black", "Edges" = "firebrick")) +
  geom_text(data = location_cor,aes(x = x_limits[1] + 0.02 * diff(x_limits),y = Inf,label = label),inherit.aes = FALSE,hjust = 0,vjust = 1.2,size = 4.5,fontface = "bold") +
  labs( x = "Date", y = "Monthly proportion", color = NULL) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom",strip.text = element_text(face = "bold"))

p_faceted

ggsave("faceted_location_observed_vs_edges_with_R.pdf",plot = p_faceted,height = 8,width = 12,units = "in")

##########


location_scatter <- location_prop_obs %>%
  inner_join(
    location_prop_edge,
    by = c("month_decimal", "location"))

# faceted scatter plot: Edges vs Donors
p_scatter <- ggplot(location_scatter, aes(x = prop_obs, y = prop_edge)) +
  geom_point(size = 2,alpha = 0.7) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed",
              color = "grey60") +
  facet_wrap(~ location, scales = "fixed") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = 0.03)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = 0.03)) +
  geom_text(data = location_cor,aes(  x = 0,  y = Inf,  label = label),
            inherit.aes = FALSE,
            hjust = 0,vjust = 1.2,size = 4.5,fontface = "bold") +
  labs(x = "Donor proportion",y = "Edge proportion") +
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold"))

# draw plot
p_scatter

# save plot
ggsave("faceted_location_edges_vs_donors_scatter.pdf",plot = p_scatter,height = 8,width = 12,units = "in")


interval_size <- 1 / 12
x_limits <- c(2021.5, 2025.33)

# donor seg propos
location_prop_obs <- plot_df %>%
  mutate(
    location = sub("\\+.*", "", location),
    location = recode(
      location,
      atlantic_flyway    = "atl",
      mississippi_flyway = "mis",
      central_flyway     = "cen",
      pacific_flyway     = "pac"
    ),
    month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, location, name = "n_obs") %>%
  group_by(month_decimal) %>%
  mutate(prop_obs = n_obs / sum(n_obs)) %>%
  ungroup() %>%
  select(month_decimal, location, prop_obs)

# edge props 

location_prop_edge <- as_tibble(beast_tree) %>%
  mutate(
    location = sub("\\+.*", "", location),
    year = 2025.33 - as.numeric(height),
    month_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, location, name = "n_edge") %>%
  group_by(month_decimal) %>%
  mutate(prop_edge = n_edge / sum(n_edge)) %>%
  ungroup() %>%
  select(month_decimal, location, prop_edge)


location_mixed <- location_prop_obs %>%
  full_join(
    location_prop_edge,
    by = c("month_decimal", "location")
  ) %>%
  replace_na(list(prop_obs = 0, prop_edge = 0)) %>%
  mutate(
    donor_bin = factor(
      ifelse(prop_obs > 0, 1, 0),
      levels = c(0, 1),
      labels = c("Absent", "Present")
    )
  )

location_stats <- location_mixed %>%
  group_by(location) %>%
  summarise(
    R = if (length(unique(donor_bin)) > 1)
      cor(as.numeric(donor_bin) - 1, prop_edge, method = "pearson")
    else NA_real_,
    p_val = if (length(unique(donor_bin)) > 1)
      t.test(prop_edge ~ donor_bin)$p.value
    else NA_real_,
    .groups = "drop") %>%
  mutate(label = paste0(
      "R = ", round(R, 3),
      "\np = ", formatC(p_val, format = "e", digits = 2)))

# vis

p_mixed_location <- ggplot(
  location_mixed,
  aes(x = prop_edge, y = donor_bin)) +
  geom_jitter(height = 0.15,width = 0,size = 6,alpha = 0.7) +
  facet_wrap(~ location, scales = "fixed") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),expand = expansion(mult = 0.03)) +
  geom_text(data = location_stats,aes(x = 0.61, y = "Present", label = label),inherit.aes = FALSE,hjust = 0,vjust = 0,size = 4.2,fontface = "bold") +
  labs(x = "Edge proportion",y = "Donor presence") +
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold"))

p_mixed_location

ggsave("faceted_location_edges_vs_donor_binary_mixed.pdf",plot = p_mixed_location,height = 8,width = 12,units = "in")




