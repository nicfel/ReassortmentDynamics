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


setwd("~/Dropbox/targeted-beast-DTA-HPAI-reassortment/MKJ-targeted/results/v2/MKJ-targeted-order-duckgoose_v2/")
most_recent_date <- 2025.33   # 2025-05-01 mrsd
tree_folder <- "MCC/"
metadata_file <- "metadata-with-clade.tsv"
# metadta file comes from moncla lab github repo for nextstrain

order_colors <- c(
  acc   = "#0072B2",
  cha   = "#009E73",
  duck  = "#56B4E9",
  gal   = "#D55E00",
  goose = "#CC79A7",
  nhm   = "#999933",
  pas   = "#E69F00",
  stri  = "#882255",
  swan  = "#117733" 
)


# Read metadata
metadata <- read_tsv(metadata_file, col_types = cols())


# Function to calculate TMRCA per genoflu group
process_tree <- function(tree_file, metadata, most_recent_date) {
  cat("Processing:", tree_file, "\n")
  tr <- read.beast(tree_file)
  
  # Extract tip labels and clean strain names
  tip_data <- tibble(label = tr@phylo$tip.label) %>%
    mutate(strain = str_split_fixed(label, "\\|", 2)[,1])
  
  # Join with metadata
  tip_data <- tip_data %>%
    left_join(metadata, by = "strain")
  
  # Filter out missing genoflu
  present <- tip_data %>% filter(!is.na(genoflu))
  
  # For each genoflu group, find MRCA node number
  # make sure MRCA is numeric (node index in phylo)
  res <- present %>%
    group_by(genoflu) %>%
    summarise(
      mrca_node = MRCA(tr@phylo, label),
      .groups = "drop") %>%
    mutate(mrca_node = as.integer(mrca_node)) %>%
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
      
      tmrca_decimal_date_low = most_recent_date - height_low,
      tmrca_decimal_date_high = most_recent_date - height_high,
      ord = ord,
      ord_prob = as.numeric(ord.prob)) %>%
    select(genoflu, tmrca_height_years_bp, tmrca_decimal_date, height_range, height_low, height_high, 
           tmrca_decimal_date_low, tmrca_decimal_date_high, ord, ord_prob)
  
  res$tree_file <- basename(tree_file)
  return(res)
}

# Process all trees in folder
tree_files <- list.files(tree_folder, pattern = "*.tree$", full.names = TRUE)

all_results <- map_dfr(tree_files, ~ process_tree(.x, metadata, most_recent_date))


write_csv(all_results, "tmrca_results.csv")
all_results



tr_p <- read.beast("MCC/ordgooseduck_MKJ_HA.mcc.tree")

tip_data_p <- tibble(label = tr_p@phylo$tip.label) %>%
  mutate(strain = str_split_fixed(label, "\\|", 2)[,1])

treedata <- left_join(tip_data_p,metadata, by = "strain" )
genoflu_over10 <- names(which(table(treedata$genoflu) > 10))



####################################
s
segment_order <- c("NS", "MP", "NA", "NP", "HA", "PA", "PB1", "PB2")

plot_df <- all_results %>%
  filter(genoflu %in% genoflu_over10) %>%
  filter(!str_detect(genoflu, "Minor")) %>%
  filter(!str_detect(genoflu, "Unseen")) %>%
  filter(!str_detect(genoflu, "divergent")) %>%
  mutate(tree_file = str_remove_all(tree_file, "ordgooseduck_MKJ_")) %>%
  mutate(tree_file = str_remove_all(tree_file, "\\.mcc\\.tree")) %>%
  mutate(ord = sub("\\+.*", "", ord)) %>%
  mutate(tree_file = factor(tree_file, levels = segment_order))

order_colors <- c(
  acc   = "#0072B2",
  cha   = "#009E73",
  duck  = "#56B4E9",
  gal   = "#D55E00",
  goose = "#CC79A7",
  nhm   = "#999933",
  pas   = "#E69F00",
  stri  = "#882255",
  swan  = "#117733"
)

p <- ggplot()
ords <- unique(plot_df$ord)

for (i in seq_along(ords)) {
  this_ord <- ords[i]
  df_ord <- plot_df %>% filter(ord == this_ord)
  
  if (i > 1) p <- p + new_scale("alpha")
  
  p <- p +
    geom_errorbarh(
      data = df_ord,
      aes(
        y = tree_file,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high,
        alpha = ord_prob),
      color = order_colors[this_ord],
      height = 0.25,
      linewidth = 1) +
    geom_point(
      data = df_ord,
      aes(
        x = tmrca_decimal_date,
        y = tree_file,
        alpha = ord_prob),
      color = order_colors[this_ord],
      size = 5) +
    scale_alpha(
      range = c(0.3, 1),
      name = paste0(this_ord, " prob"))
  
  df_grey <- df_ord %>% filter(ord_prob < 0.5)
  
  p <- p +
    geom_errorbarh(
      data = df_grey,
      aes(
        y = tree_file,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high),
      color = "grey70",
      height = 0.25,
      linewidth = 1) +
    geom_point(data = df_grey, aes(
        x = tmrca_decimal_date,
        y = tree_file),
      color = "grey70",
      size = 5)
}

p <- p +
  facet_wrap(~ genoflu, ncol = 4) +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = "Date",
    y = "Segment") +
  theme_minimal(base_size = 14)

p

ggsave("faceted_genotype_tmcra_order.pdf",plot = p,height = 12,width = 11,units = "in")

##############################

segment_order <- c("NS", "MP", "NA", "NP", "HA", "PA", "PB1", "PB2")

plot_df <- all_results %>%
  filter(genoflu %in% genoflu_over10) %>%
  filter(!str_detect(genoflu, "Minor")) %>%
  filter(!str_detect(genoflu, "Unseen")) %>%
  filter(!str_detect(genoflu, "divergent")) %>%
  mutate(
    tree_file = str_remove_all(tree_file, "ordgooseduck_MKJ_"),
    tree_file = str_remove_all(tree_file, "\\.mcc\\.tree"),
    ord = sub("\\+.*", "", ord),
    tree_file = factor(tree_file, levels = segment_order),
    genoflu = factor(genoflu, levels = sort(unique(genoflu))),
    geno_segment = paste(genoflu, tree_file, sep = " | "))

order_colors <- c(
  acc   = "#0072B2",
  cha   = "#009E73",
  duck  = "#56B4E9",
  gal   = "#D55E00",
  goose = "#CC79A7",
  nhm   = "#999933",
  pas   = "#E69F00",
  stri  = "#882255",
  swan  = "#117733"
)

p <- ggplot()
ords <- unique(plot_df$ord)

for (i in seq_along(ords)) {
  this_ord <- ords[i]
  df_ord <- plot_df %>% filter(ord == this_ord)
  
  if (i > 1) p <- p + new_scale("alpha")
  
  p <- p +
    geom_errorbarh(
      data = df_ord,
      aes(
        y = geno_segment,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high,
        alpha = ord_prob),
      color = order_colors[this_ord],
      height = 0.2,
      linewidth = 1.2) +
    geom_point(data = df_ord,aes(
        x = tmrca_decimal_date,
        y = geno_segment,
        alpha = ord_prob),
      color = order_colors[this_ord],
      size = 2) +
    scale_alpha(
      range = c(0.3, 1),
      name = paste0(this_ord, " prob"))
  
  df_grey <- df_ord %>% filter(ord_prob < 0.5)
  
  p <- p +
    geom_errorbarh(
      data = df_grey,
      aes(
        y = geno_segment,
        xmin = tmrca_decimal_date_low,
        xmax = tmrca_decimal_date_high),
      color = "grey70",
      height = 0.2,
      linewidth = 1.2) +
    geom_point(data = df_grey,aes(x = tmrca_decimal_date, y = geno_segment), color = "grey70", size = 2)
}

p_final <- p +
  labs(
    title = "TMRCA estimates across all genotypes (excluding Minor)",
    x = "Date",
    y = NULL) +
  facet_grid(
    genoflu ~ .,
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

p_final

ggsave("figure_supp_allgeno_order.pdf",plot = p_final,height = 11,width = 11,units = "in")




##############################
# edge colorings based on the state
# Based on reweards  visualization
# use the MCC tree with edges annotated to the branches

beast_tree <- read.beast("MCC/ordgooseduck_MKJ_HA.mcc.tree")
data<-as_tibble(beast_tree)
data <- data %>% 
  mutate(ord = sub("\\+.*", "", ord))
## filter data
df<- data %>% select(ord,height)

# input date of earliest seq
df <- df %>% mutate(year = 2025.33 - as.numeric(height)) %>%
  mutate(ord = sub("\\+.*", "", ord))
## get the propotrion of rewards

interval_size <- 1/52   # weekly decimal interval

df_long_weekly <- df %>%
  mutate(
    week_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  group_by(week_decimal, ord) %>%
  summarise(total_freq = n(), .groups = "drop")

all_weeks <- sort(unique(df_long_weekly$week_decimal))
all_ords  <- sort(unique(df_long_weekly$ord))

df_long_weekly_complete <- df_long_weekly %>%
  complete(week_decimal = all_weeks, ord = all_ords, fill = list(total_freq = 0))

df_long_weekly_percent <- df_long_weekly_complete %>%
  group_by(week_decimal) %>%
  mutate(
    prop = if (sum(total_freq) == 0) 0 else total_freq / sum(total_freq)) %>%
  ungroup()

p_week_decimal <- ggplot(df_long_weekly_percent,
                         aes(x = week_decimal, y = prop, fill = ord)) +
  geom_area(size = 0.2, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0,0)) +
  xlim(2022, 2025.33) +
  labs(
    title = "% Edges by Host Order (Weekly)",
    x = "Year",
    y = "Proportion",
    fill = "Order") +
  theme_minimal(base_size = 18) +
  scale_fill_manual(values = order_colors)

p_week_decimal

ggsave("stacked_week_order.pdf", plot = p_week_decimal, height = 8.5, width = 12, units = "in")



###########
# plot with the annoations for tmcra of the genotyope

p_week_decimal2 <- p_week_decimal +
  geom_vline(
    data = plot_df,
    aes(xintercept = tmrca_decimal_date, color = ord),
    linetype = "dashed",
    linewidth = 0.7,
    alpha = 0.8) +
  geom_text(
    data = plot_df,
    aes(x = tmrca_decimal_date,label = genoflu,y= 1,color = ord),
    angle = 90,
    vjust = -0.2,
    size = 4,
    fontface = "bold") +
  scale_color_manual(values = order_colors, name = "Order") +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)) +
  coord_cartesian(clip = "off") 


p_week_decimal2

# Save the plot
ggsave("stacked_week_order-withgenos.pdf", plot = p_week_decimal2, height = 8.5, width = 12, units = "in")



###########################
# plot of genotypes over time vs the proportion of orders
# filter out those which have posterior probabilty lower than 0.5

plot_df <- all_results %>%
  filter(genoflu %in% genoflu_over10) %>%
  filter(!str_detect(genoflu, "Minor")) %>%
  filter(!str_detect(genoflu, "Unseen")) %>%
  filter(!str_detect(genoflu, "divergent")) %>%
  mutate(
    tree_file = str_remove_all(tree_file, "ordgooseduck_MKJ_"),
    tree_file = str_remove_all(tree_file, ".mcc.tree"),
    ord = sub("\\+.*", "", ord),
    tree_file = factor(tree_file, levels = segment_order),
    genoflu = factor(genoflu, levels = sort(unique(genoflu))),
    geno_segment = paste(genoflu, tree_file, sep = " | ")
  )


order_colors <- c(
  acc   = "#0072B2",
  cha   = "#009E73",
  duck  = "#56B4E9",
  gal   = "#D55E00",
  goose = "#CC79A7",
  nhm   = "#999933",
  pas   = "#E69F00",
  stri  = "#882255",
  swan  = "#117733" 
)

plot_df <- plot_df %>%
  mutate(date = date_decimal(tmrca_decimal_date))

# Count occurrences per month per 'ord'
monthly_counts <- plot_df %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(ord, year_month) %>%
  summarise(count = n(), .groups = "drop")

# Complete time series for each 'ord'
all_months <- seq(min(monthly_counts$year_month), max(monthly_counts$year_month), by = "month")
ords <- unique(monthly_counts$ord)
full_ts <- expand.grid(ord = ords, year_month = all_months)
monthly_counts <- full_ts %>%
  left_join(monthly_counts, by = c("ord", "year_month")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  arrange(ord, year_month)

monthly_counts <- monthly_counts %>%
  group_by(ord) %>%
  mutate(rolling_avg = rollmean(count, k = 4, fill = NA, align = "right")) %>%
  ungroup()

monthly_counts$year_month <- as.Date(monthly_counts$year_month)

geno_monthly <- plot_df %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(year_month) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  arrange(year_month) %>%
  mutate(rolling_total = rollmean(total_count, k = 3, fill = NA, align = "right"))

geno_monthly$year_month <- as.Date(geno_monthly$year_month)
geno_monthly <- geno_monthly %>%
  mutate(Genotype = "Genotype segments")

# Plot
sp <- ggplot(monthly_counts, aes(x = year_month, y = rolling_avg, color = ord)) +
  geom_line(size = 2) +
  geom_line(data = geno_monthly, aes(x = year_month, y = rolling_total, linetype = Genotype, color = Genotype), 
            size = 1) +
  scale_color_manual(values = order_colors) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b\n%Y") +
  scale_linetype_manual(
    values = c("Genotype segments" = "dashed")) +
  labs(
    title = "Number of segments from Host ",
    x = "Date",
    y = "Rolling 3-Month Average Count",
    color = "Host") +
  theme_minimal(base_size = 18)


sp

ggsave("segemnts_genotype_host.pdf", plot = sp, height = 8.5, width = 12, units = "in")

