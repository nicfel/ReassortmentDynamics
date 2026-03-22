library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(progressr)
library(zoo)
# script to do observed/expected ratio for orders and flyawys states 

set.seed(1738)

n_draws <- 10000
n_boot  <- 100

handlers(global = TRUE)
handlers("txtprogressbar")

interval_size <- 1 / 12
x_limits <- c(2021.5, 2025)
area_scale <- 1
two_week_window <- 14 / 365

##########
# O/E test with twsided p value 

compute_oe_pvals <- function(df, group_var, draw_var, n_draws = 1000) {
  
  observed_counts <- df %>%
    group_by(.data[[group_var]]) %>%
    summarise(observed_n = sum(observed_n), .groups = "drop")
  
  sim_results <- map_dfr(1:n_draws, function(sim_id) {
    df %>%
      mutate(
        sim = sim_id,
        sim_value = map_chr(edge_props, function(ep) {
          sample(ep[[draw_var]], size = 1, replace = TRUE, prob = ep$prop)
        })
      )
  })
  
  sim_counts <- sim_results %>%
    count(sim, value = sim_value, name = "sim_n") %>%
    rename(!!group_var := value) %>%
    filter(.data[[group_var]] %in% observed_counts[[group_var]])
  
  expected_summary <- sim_counts %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      expected_mean  = mean(sim_n),
      expected_lower = quantile(sim_n, 0.025),
      expected_upper = quantile(sim_n, 0.975),
      .groups = "drop"
    )
  
  # Two-sided p-values
  pvals <- sim_counts %>%
    left_join(observed_counts, by = group_var) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      p_upper = (sum(sim_n >= observed_n) + 1) / (n_draws + 1),
      p_lower = (sum(sim_n <= observed_n) + 1) / (n_draws + 1),
      p_value = pmin(1, 2 * pmin(p_upper, p_lower)),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = formatC(p_value, format = "e", digits = 2),
      signif  = ifelse(p_value < 0.05, "*", "")
    )
  
  observed_counts %>%
    left_join(expected_summary, by = group_var) %>%
    left_join(pvals, by = group_var) %>%
    mutate(
      oe_ratio = observed_n / expected_mean,
      oe_lower = observed_n / expected_upper,
      oe_upper = observed_n / expected_lower
    )
}


# boostrap function 

bootstrap_oe <- function(df, group_var, draw_var, n_boot, n_draws) {
  
  df <- df %>%
    mutate(genotype = sub(" \\|.*", "", geno_segment))
  
  genotypes <- unique(df$genotype)
  N_geno <- length(genotypes)
  
  with_progress({
    
    p <- progressor(steps = n_boot)
    
    boot_results <- map_dfr(1:n_boot, function(b) {
      
      p(sprintf("Bootstrap %d", b))
      
      sampled_genotypes <- sample(genotypes, size = N_geno, replace = TRUE)
      
      boot_sample <- map_dfr(sampled_genotypes, function(g) {
        df %>% filter(genotype == g)
      })
      
      res <- compute_oe_pvals(
        boot_sample,
        group_var = group_var,
        draw_var  = draw_var,
        n_draws   = n_draws
      )
      
      res$boot_id <- b
      res
    })
  })
  
  boot_results %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      oe_lower_boot = quantile(oe_ratio, 0.025),
      oe_upper_boot = quantile(oe_ratio, 0.975),
      .groups = "drop"
    )
}


# Get segment edge proportions

build_segment_edge_props <- function(plot_df, edge_table, group_var) {
  
  segments_interval <- plot_df %>%
    select(geno_segment, !!sym(group_var), tmrca_decimal_date) %>%
    mutate(
      interval_low  = tmrca_decimal_date - two_week_window / 2,
      interval_high = tmrca_decimal_date + two_week_window / 2
    )
  
  segments_interval %>%
    rowwise() %>%
    mutate(
      edges_in_interval = list(
        edge_table %>%
          filter(edge_decimal_date >= interval_low,
                 edge_decimal_date <= interval_high)
      ),
      observed_n = ifelse(
        nrow(edges_in_interval) > 0 &
          any(edges_in_interval[[group_var]] == .data[[group_var]]),
        1, 0
      )
    ) %>%
    ungroup() %>%
    filter(lengths(edges_in_interval) > 0) %>%
    rowwise() %>%
    mutate(
      edge_props = list({
        edges <- edges_in_interval
        if (nrow(edges) == 0) NULL else
          edges %>%
          count(!!sym(group_var), name = "n") %>%
          mutate(prop = n / sum(n))
      })
    ) %>%
    ungroup() %>%
    filter(!map_lgl(edge_props, is.null))
}


# visaualize plot

plot_panel_C <- function(df, group_var, fill_colors, title, xlab) {
  
  df %>%
    mutate(group = .data[[group_var]],
           group = factor(group, levels = group[order(oe_ratio)])) %>%
    ggplot(aes(x = group, y = oe_ratio, fill = group)) +
    geom_errorbar(aes(ymin = oe_lower_boot, ymax = oe_upper_boot), width = 0.25, linewidth = 0.8) +
    geom_point(size = 7, shape = 21, color = "black", stroke = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8) +
    geom_text(aes(label = paste0(p_label, signif), y = oe_upper_boot * 1.05), size = 4) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title,x = xlab,y = "Observed / Expected") +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(face = "bold", size = 20),legend.position = "none")
}


############################################
############################################
# ORDERS ANALYSIS 

setwd("~/Dropbox/targeted-beast-DTA-HPAI-reassortment/MKJ-targeted/results/v2/MKJ-targeted-order-duckgoose_v2/")
plot_df <-  read.csv("plot_df_order.csv")
beast_tree <- read.beast("MCC/ordgooseduck_MKJ_HA.mcc.tree")

# EDGE TABLE
edge_table_ord <- as_tibble(beast_tree) %>%
  mutate(
    ord = sub("\\+.*", "", ord),
    edge_decimal_date = 2025.33 - as.numeric(height)) %>%
  select(ord, edge_decimal_date)

# SEGMENT TABLE
segment_edge_props <- build_segment_edge_props(
  plot_df,
  edge_table_ord,
  "ord"
)

# OE + p-values
orders_results <- compute_oe_pvals(segment_edge_props,"ord","ord",n_draws)

# Bootstrap CI
orders_boot <- bootstrap_oe(segment_edge_props,"ord","ord",n_boot,n_draws)

final_orders <- orders_results %>%
  left_join(orders_boot, by = "ord")


# Visualization orders analysis 

p_orders_C <- plot_panel_C( final_orders, "ord", order_colors, "C) Orders (Bootstrap CI + Two-Sided Test)", "Host order")


############################################
############################################
## FLYWAY ANALYSIS
setwd("~/Dropbox/targeted-beast-DTA-HPAI-reassortment/h5-data-updates-main-reassortmentprev/h5nx/")
plot_df <- read.csv("plot_df_flyway.csv")
beast_tree <- read.beast("MCC/HPAI_HA_northamerica_targeted_dta.mcc.tree")

# EDGE TABLE
edge_table_fly <- as_tibble(beast_tree) %>%
  mutate(
    location = sub("\\+.*", "", location),
    edge_decimal_date = 2025.33 - as.numeric(height)
  ) %>%
  select(location, edge_decimal_date)

edge_table_fly$location <- recode(
  edge_table_fly$location,
  pac = "pacific_flyway",
  mis = "mississippi_flyway",
  cen = "central_flyway",
  atl = "atlantic_flyway"
)

# SEGMENT TABLE
segment_edge_props_flyway <- build_segment_edge_props(plot_df,edge_table_fly,"location")

# OE + p-values
flyway_results <- compute_oe_pvals(segment_edge_props_flyway,"location","location",n_draws)

# Bootstrap CI
flyway_boot <- bootstrap_oe(segment_edge_props_flyway,"location","location",n_boot,n_draws)

final_flyway <- flyway_results %>%
  left_join(flyway_boot, by = "location")


p_flyway_C <- plot_panel_C(final_flyway,"location",flyway_colors,"C) Flyways (Bootstrap CI + Two-Sided Test)","Flyway")


write.csv(final_orders,  "final_results_orders.csv",  row.names = FALSE)
write.csv(final_flyway, "final_results_flyway.csv", row.names = FALSE)

ggsave("panel_C_orders.pdf",  p_orders_C,  width = 8, height = 6)
ggsave("panel_C_flyways.pdf", p_flyway_C, width = 8, height = 6)