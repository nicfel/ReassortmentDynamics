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
      aes( y = tree_file, xmin = tmrca_decimal_date_low, xmax = tmrca_decimal_date_high, alpha = ord_prob),
      color = order_colors[this_ord],
      height = 0.25,
      linewidth = 1) +
    geom_point(data = df_ord,aes(x = tmrca_decimal_date,y = tree_file,alpha = ord_prob),color = order_colors[this_ord],size = 5) +
    scale_alpha(range = c(0.3, 1),name = paste0(this_ord, " prob"))
  
  df_grey <- df_ord %>% filter(ord_prob < 0.5)
  
  p <- p +
    geom_errorbarh(
      data = df_grey,
      aes(y = tree_file,xmin = tmrca_decimal_date_low,xmax = tmrca_decimal_date_high),
      color = "grey70",
      height = 0.25,
      linewidth = 1) +
    geom_point(data = df_grey, aes(x = tmrca_decimal_date,y = tree_file),color = "grey70",size = 5)
}

p <- p +
  facet_wrap(~ genoflu, ncol = 4) +
  labs( title = "TMRCA estimates across all genotypes (excluding Minor)", x = "Date", y = "Segment") +
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

write.csv(plot_df,"plot_df_order.csv", row.names = FALSE)

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
      aes(y = geno_segment,xmin = tmrca_decimal_date_low,xmax = tmrca_decimal_date_high,alpha = ord_prob),
      color = order_colors[this_ord],
      height = 0.2,
      linewidth = 1.2) +
    geom_point(data = df_ord,aes(x = tmrca_decimal_date,y = geno_segment,alpha = ord_prob),
      color = order_colors[this_ord],
      size = 2) +
    scale_alpha(range = c(0.3, 1),name = paste0(this_ord, " prob"))
  
  df_grey <- df_ord %>% filter(ord_prob < 0.5)
  
  p <- p +
    geom_errorbarh(data = df_grey,aes(y = geno_segment,xmin = tmrca_decimal_date_low,xmax = tmrca_decimal_date_high),color = "grey70",height = 0.2,linewidth = 1.2) +
    geom_point(data = df_grey,aes(x = tmrca_decimal_date, y = geno_segment), color = "grey70", size = 2)
}

p_final <- p +
  labs(title = "TMRCA estimates across all genotypes (excluding Minor)",x = "Date",y = NULL) +
  facet_grid( genoflu ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal(base_size = 14) +
  theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"), strip.background = element_blank(), panel.spacing.y = unit(0.025, "lines")) +
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
  mutate(week_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
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

p_week_decimal <- ggplot(df_long_weekly_percent, aes(x = week_decimal, y = prop, fill = ord)) +
  geom_area(size = 0.2, alpha = 0.9) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = c(0,0)) +
  xlim(2021.5, 2025.33) +
  labs( title = "% Edges by Host Order (Weekly)", x = "Year", y = "Proportion", fill = "Order") +
  theme_minimal(base_size = 18) +
  scale_fill_manual(values = order_colors)

p_week_decimal

ggsave("stacked_week_order.pdf", plot = p_week_decimal, height = 8.5, width = 12, units = "in")



###########
# plot with the annoations for tmcra of the genotyope

p_week_decimal2 <- p_week_decimal +
  geom_vline(data = plot_df,aes(xintercept = tmrca_decimal_date, color = ord),linetype = "dashed",linewidth = 0.7,alpha = 0.8) +
  geom_text(data = plot_df,aes(x = tmrca_decimal_date,label = genoflu,y= 1,color = ord),angle = 90,vjust = -0.2,size = 4,fontface = "bold") +
  scale_color_manual(values = order_colors, name = "Order") +
  guides(fill = guide_legend(order = 1),color = guide_legend(order = 2)) +
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
  geom_line(data = geno_monthly, aes(x = year_month, y = rolling_total, linetype = Genotype, color = Genotype), size = 1) +
  scale_color_manual(values = order_colors) +
  scale_x_date(date_breaks = "4 month", date_labels = "%b\n%Y") +
  scale_linetype_manual(values = c("Genotype segments" = "dashed")) +
  labs(title = "Number of segments from Host ",x = "Date",y = "Rolling 3-Month Average Count",color = "Host") +
  theme_minimal(base_size = 18)


sp

ggsave("segemnts_genotype_host.pdf", plot = sp, height = 8.5, width = 12, units = "in")


#################################
interval_size <- 1 / 12 
x_limits <- c(2021.5, 2025)
area_scale <- 1


geno_monthly <- plot_df %>%
  mutate(month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, name = "geno_n")

# panel A - donors of segements

ord_prop <- plot_df %>%
  mutate(
    ord = sub("\\+.*", "", ord),
    month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_ord") %>%
  complete(month_decimal = sort(unique(month_decimal)),ord = sort(unique(ord)),fill = list(n_ord = 0)) %>%
  group_by(month_decimal) %>%
  mutate(prop = if (sum(n_ord) == 0) 0 else n_ord / sum(n_ord)) 

ord_scaled <- ord_prop %>%
  left_join(geno_monthly, by = "month_decimal") %>%
  mutate(geno_n = replace_na(geno_n, 0),scaled_height = prop * geno_n * area_scale)

p_panel_A <- ggplot() +
  geom_area(data = ord_scaled,aes(x = month_decimal, y = scaled_height, fill = ord),alpha = 0.9,linewidth = 0.2) +
  geom_line(data = geno_monthly,aes(x = month_decimal, y = geno_n),linewidth = 1.3,color = "black") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = order_colors) +
  labs(title = "A) Donor segment proportion",y = "Genotype segements",x = NULL,fill = "Order") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20))

p_panel_A
# panel b - edges of the phylogeny 

edge_prop <- as_tibble(beast_tree) %>%
  mutate(ord = sub("\\+.*", "", ord),year = 2025.33 - as.numeric(height),month_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_edges")



all_months <- seq(min(edge_prop$month_decimal), max(edge_prop$month_decimal), by = interval_size)
all_ords <- sort(unique(edge_prop$ord))

edge_prop_complete <- edge_prop %>%
  complete(month_decimal = all_months,ord = all_ords,fill = list(n_edges = 0)) %>%
  group_by(month_decimal) %>%
  mutate(prop = if (sum(n_edges) == 0) 0 else n_edges / sum(n_edges)) %>%
  ungroup()

edge_scaled <- edge_prop_complete %>%
  left_join(geno_monthly, by = "month_decimal") %>%
  arrange(month_decimal) %>%
  group_by(ord) %>%
  mutate(geno_n = na.approx(geno_n, x = month_decimal, na.rm = FALSE, rule = 2), scaled_height = prop * geno_n * area_scale) %>%
  ungroup()

p_panel_B <- ggplot() +
  geom_area( data = edge_scaled, aes(x = month_decimal, y = scaled_height, fill = ord), alpha = 0.9, linewidth = 0.2) +
  geom_line(data = geno_monthly,aes(x = month_decimal, y = geno_n),linewidth = 1.3,color = "black") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = order_colors) +
  labs( title = "B) Edges proportion", y = "Genotype segments", x = "Year", fill = "Order") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(face = "bold", size = 20))


p_panel_B


final_figure <- p_panel_A / p_panel_B +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

final_figure

ggsave( "geno_edges_vs_donors_AB.pdf", plot = final_figure, height = 12, width = 12, units = "in")



###########
# correlation calculation 

interval_size <- 1 / 12
x_limits <- c(2021.5, 2025.33)
cor_method <- "spearman"

# donor segments calc
ord_prop_obs <- plot_df %>%
  mutate(
    ord = sub("\\+.*", "", ord),
    month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_obs") %>%
  group_by(month_decimal) %>%
  mutate(prop_obs = n_obs / sum(n_obs)) %>%
  ungroup() %>%
  select(month_decimal, ord, prop_obs)


# edge proportions
ord_prop_edge <- as_tibble(beast_tree) %>%
  mutate(
    ord = sub("\\+.*", "", ord),
    year = 2025.33 - as.numeric(height),
    month_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_edge") %>%
  group_by(month_decimal) %>%
  mutate(prop_edge = n_edge / sum(n_edge)) %>%
  ungroup() %>%
  select(month_decimal, ord, prop_edge)


# format data
ord_long <- ord_prop_obs %>%
  full_join(ord_prop_edge,
            by = c("month_decimal", "ord")) %>%
  replace_na(list(prop_obs = 0, prop_edge = 0)) %>%
  pivot_longer(
    cols = c(prop_obs, prop_edge),
    names_to = "source",
    values_to = "prop") %>%
  mutate(
    source = recode(source, prop_obs = "Donor", prop_edge = "Edges"))


# calc spearman corrl
ord_cor <- ord_prop_obs %>%
  inner_join(ord_prop_edge,
             by = c("month_decimal", "ord")) %>%
  group_by(ord) %>%
  summarise(
    R = if (n() >= 1)
      cor(prop_obs, prop_edge, method = cor_method), .groups = "drop") %>%
  mutate(label = paste0("R = ", round(R, 2)))


p_faceted <- ggplot(ord_long,aes(x = month_decimal, y = prop, color = source)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  facet_wrap(~ ord, scales = "fixed") +
  scale_x_continuous(limits = x_limits,expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),expand = c(0, 0)) +
  scale_color_manual(values = c("Donor" = "black", "Edges" = "firebrick")) +
  geom_text(data = ord_cor, aes(x = x_limits[1] + 0.02 * diff(x_limits),y = Inf,label = label), inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 4.5, fontface = "bold") +
  labs( x = "Date", y = "Monthly proportion", color = NULL) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom",strip.text = element_text(face = "bold"))

p_faceted

ggsave("faceted_ord_observed_vs_edges_with_R.pdf",plot = p_faceted,height = 8,width = 12,units = "in")

##########


ord_scatter <- ord_prop_obs %>%
  inner_join(
    ord_prop_edge,
    by = c("month_decimal", "ord"))

# faceted scatter plot: Edges vs Donors
p_scatter <- ggplot(ord_scatter, aes(x = prop_obs, y = prop_edge)) +
  geom_point(size = 2,alpha = 0.7) +
  geom_abline(slope = 1,intercept = 0,linetype = "dashed",color = "grey60") +
  facet_wrap(~ ord, scales = "fixed") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = 0.03)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = 0.03)) +
  geom_text(data = ord_cor,aes(  x = 0,  y = Inf,  label = label), inherit.aes = FALSE, hjust = 0,vjust = 1.2,size = 4.5,fontface = "bold") +
  labs(x = "Donor proportion",y = "Edge proportion") +
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold"))

p_scatter

ggsave("faceted_ord_edges_vs_donors_scatter.pdf",plot = p_scatter,height = 8,width = 12,units = "in")




############################
# binary correlations 

interval_size <- 1 / 12
x_limits <- c(2021.5, 2025.33)

# donor segs
ord_prop_obs <- plot_df %>%
  mutate(
    ord = sub("\\+.*", "", ord),
    month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_obs") %>%
  group_by(month_decimal) %>%
  mutate(prop_obs = n_obs / sum(n_obs)) %>%
  ungroup() %>%
  select(month_decimal, ord, prop_obs)

# edge props 
ord_prop_edge <- as_tibble(beast_tree) %>%
  mutate(
    ord = sub("\\+.*", "", ord),
    year = 2025.33 - as.numeric(height),
    month_decimal = (floor(year / interval_size) + 0.5) * interval_size
  ) %>%
  count(month_decimal, ord, name = "n_edge") %>%
  group_by(month_decimal) %>%
  mutate(prop_edge = n_edge / sum(n_edge)) %>%
  ungroup() %>%
  select(month_decimal, ord, prop_edge)


ord_binary <- ord_prop_obs %>%
  full_join(
    ord_prop_edge,
    by = c("month_decimal", "ord")
  ) %>%
  replace_na(list(prop_obs = 0, prop_edge = 0)) %>%
  mutate(
    bin_obs  = as.integer(prop_obs  > 0),
    bin_edge = as.integer(prop_edge > 0)
  )


ord_stats <- ord_binary %>%
  group_by(ord) %>%
  summarise(
    R = if (n() >= 2)
      cor(bin_obs, bin_edge, method = "pearson")
    else NA_real_,
    fisher_p = {
      tab <- table(bin_obs, bin_edge)
      if (all(dim(tab) == c(2, 2)))
        fisher.test(tab)$p.value
      else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "R = ", round(R, 2),
      "\nFisher p = ", formatC(fisher_p, format = "e", digits = 2)
    )
  )


p_scatter_bin <- ggplot(ord_binary,aes(x = bin_obs, y = bin_edge)) +
  geom_jitter(width = 0.08,height = 0.08,size = 2,alpha = 0.7) +
  facet_wrap(~ ord, scales = "fixed") +
  scale_x_continuous(breaks = c(0, 1),labels = c("Absent", "Present"),limits = c(-0.2, 1.2)) +
  scale_y_continuous(breaks = c(0, 1),labels = c("Absent", "Present"),limits = c(-0.2, 1.2)) +
  geom_text(data = ord_stats,aes(x = -0.15, y = -0.15, label = label),inherit.aes = FALSE,hjust = 0,vjust = 0,size = 4.2,fontface = "bold") +
  labs(x = "Donor present",y = "Edge present") +
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold"))

p_scatter_bin

ggsave("faceted_ord_binary_scatter_with_fisher.pdf",plot = p_scatter_bin,height = 8,width = 12,units = "in")



# settings
interval_size <- 1 / 12
x_limits <- c(0, 1)

# donor segs
ord_prop_obs <- plot_df %>%
  mutate(ord = sub("\\+.*", "", ord),month_decimal = (floor(tmrca_decimal_date / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_obs") %>%
  group_by(month_decimal) %>%
  mutate(prop_obs = n_obs / sum(n_obs)) %>%
  ungroup() %>%
  select(month_decimal, ord, prop_obs)

# edge props
ord_prop_edge <- as_tibble(beast_tree) %>%
  mutate( ord = sub("\\+.*", "", ord), year = 2025.33 - as.numeric(height), month_decimal = (floor(year / interval_size) + 0.5) * interval_size) %>%
  count(month_decimal, ord, name = "n_edge") %>%
  group_by(month_decimal) %>%
  mutate(prop_edge = n_edge / sum(n_edge)) %>%
  ungroup() %>%
  select(month_decimal, ord, prop_edge)

###
ord_mixed <- ord_prop_obs %>%
  full_join(
    ord_prop_edge,
    by = c("month_decimal", "ord")
  ) %>%
  replace_na(list(prop_obs = 0, prop_edge = 0)) %>%
  mutate(
    donor_bin = factor(
      ifelse(prop_obs > 0, 1, 0),
      levels = c(0, 1),
      labels = c("Absent", "Present")
    )
  )


ord_stats <- ord_mixed %>%
  group_by(ord) %>%
  summarise(
    R = if (length(unique(donor_bin)) > 1)
      cor(as.numeric(donor_bin) - 1, prop_edge, method = "pearson")
    else NA_real_,
    p_val = if (length(unique(donor_bin)) > 1)
      t.test(prop_edge ~ donor_bin)$p.value
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "R = ", round(R, 2),
      "\np = ", formatC(p_val, format = "e", digits = 2)
    )
  )

# Vis
p_mixed <- ggplot(ord_mixed,aes(x = prop_edge, y = donor_bin)) +
  geom_jitter(height = 0.15,width = 0,size = 2,alpha = 0.7) +
  facet_wrap(~ ord, scales = "fixed") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),limits = x_limits,expand = expansion(mult = 0.03)) +
  geom_text(data = ord_stats,aes(x = 0.55, y = "Present", label = label),inherit.aes = FALSE,hjust = 0,vjust = 0,size = 4.2,fontface = "bold") +
  labs(x = "Edge proportion",y = "Donor presence") +
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold"))

p_mixed

ggsave("faceted_ord_mixed_edges_vs_donor_binary.pdf",plot = p_mixed,height = 8,width = 12,units = "in")




