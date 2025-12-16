library(ggplot2)
library(dplyr)
library(scales)
library(zoo)
library(patchwork)
library(ggplot2)
library(tidyr)
library(readr)
library(slider)

# sequences used in this analysis 
setwd("~/Dropbox/targeted-beast-DTA-HPAI-reassortment/h5-data-updates-main-reassortmentprev/h5nx/")
meta <- read.csv("H5N1-2021-2025_northamerica_hostorders.csv")

# Format detections of HPAI H5 in wild, domestic, mammal, cattle and combine
# Output number of detections total over time as well as by host type
# download datasets from USDA APHIS: (https://www.aphis.usda.gov/aphis/ourfocus/animalhealth/animal-disease-information/avian/avian-influenza/)
# Latest update:  2025-06-16
# for commercial flocks dataframe remove first row with "Control Area released" and collapse empty cells for columns

setwd("code/metadata/")

#### read in dataframes and format dates, make week column for CDC week

wildbirds <- read.csv("raw_detection_data_update20260616//wild_bird_detections.csv", header = TRUE)
wildbirds <- wildbirds %>%
  mutate(Collection.Date = ifelse(Collection.Date == "null", NA, Collection.Date)) %>%
  mutate(Collection.Date = ifelse(Collection.Date == "Unknown", NA, Collection.Date)) %>%
  mutate(date_det = mdy(Collection.Date)) %>%
  mutate(date_det2 = mdy(Date.Detected)) %>%
  mutate(date_det = coalesce(date_det, date_det2)) %>%
  mutate(week = floor_date(date_det, "week",  week_start = 7))



########
# add in the Canadian detection data: 
# Downloaded on 2025-12-11.
canadian_wildbird <- read.csv("raw_detection_data_update20260616/canadian_wildbird_det.csv", header = TRUE)

canadian_wildbird <- canadian_wildbird %>%
  mutate(Collection.Date = ifelse(Collection.Date == "null", NA, Collection.Date)) %>%
  mutate(Collection.Date = ifelse(Collection.Date == "Unknown", NA, Collection.Date)) %>%
  mutate(date_det = mdy(Collection.Date)) %>%
  mutate(week = floor_date(date_det, "week",  week_start = 7)) %>%
  mutate(Bird.Species = Common.Name) %>%
  mutate(State = Province)


wildbirds <- bind_rows(
  select(wildbirds, State, Bird.Species, date_det,week),
  select(canadian_wildbird, State, Bird.Species, date_det,week)
)

speciesfixs <- read.csv("species.csv")
speciesfixs <- speciesfixs %>% distinct(common_name_correction, .keep_all = TRUE)

wildbirds$Bird.Species <- gsub(" ", "_", tolower(wildbirds$Bird.Species))
wildbirds$Bird.Species <- gsub("'", "", wildbirds$Bird.Species)  # Remove single quotes
wildbirds$Bird.Species <- gsub("_\\(unidentified\\)", "", wildbirds$Bird.Species)  # Remove the string "_(unidentified)"
wildbirds <- merge(wildbirds,speciesfixs, by.x ="Bird.Species", by.y ="common_name_correction", all.x = TRUE)

species <-  read.csv("order_condensed.csv")
wildbirds <- merge(wildbirds,species, by ="order", all.x = TRUE)

wildbirds <- wildbirds %>%
  mutate(order_dg = case_when(
    broad %in% c("duck", "goose") ~ broad,
    TRUE ~ order_condensed
  ))


mammals <- read.csv("raw_detection_data_update20260616//mammal_detections.csv", header = TRUE)
mammals <- mammals %>%
  mutate(date_det = mdy(date_collected)) %>%
  mutate(week = floor_date(date_det, "week",  week_start = 7))

# for commercial flocks + backyard bird dataframe remove first row with "Control Area released" and save as UTF-8 CSV format, for some reason it has embedded null values 
dombirds <- read.csv("raw_detection_data_update20260616//domestic_bird_detections.csv", header = TRUE)
dombirds <- dombirds  %>%
  mutate(date_det = as.Date(Confirmed, format = "%d-%b-%y")) %>%
  select(date_det, State, County.Name, Special.Id, Production) %>%
  mutate(week = floor_date(date_det, "week",  week_start = 7))


# open sheet in excel and save as UTF-8 CSV
cattle <- read.csv("raw_detection_data_update20260616//cattle_detections.csv", header = TRUE)
cattle <- cattle %>%
  mutate(date_det = mdy(Confirmed)) %>%
  mutate(week = floor_date(date_det, "week",  week_start = 7))

#########
# Wild bird counts 
wildbird_count_date <- wildbirds %>%
  group_by(date_det, order_dg) %>%
  summarise(wildbird_count = n(), .groups = "drop")

wildbird_count_date_wide <- wildbird_count_date %>%
  tidyr::pivot_wider(
    names_from = order_dg,
    values_from = wildbird_count,
    values_fill = 0
  )


wildbird_count_week <- wildbirds %>%
  group_by(week, order_dg) %>%
  summarise(wildbird_count = n(), .groups = "drop")


wildbird_count_week_wide <- wildbird_count_week %>%
  tidyr::pivot_wider(
    names_from = order_dg,
    values_from = wildbird_count,
    values_fill = 0)


# domestic/commerical flock data counts
domes_count_date <- dombirds %>%
  group_by(date_det) %>%
  summarize(domestic_count = n())

domes_count_week <- dombirds %>%
  group_by(week) %>%
  summarize(domestic_count = n())

# mammal data counts
mammals_count_date <- mammals %>%
  group_by(date_det) %>%
  summarize(mammal_count = n())

mammals_count_week <- mammals %>%
  group_by(week) %>%
  summarize(mammal_count = n())

# cattle data coutns

cattle_count_date <- cattle %>%
  group_by(date_det) %>%
  summarize(cattle_count = n())

cattle_count_week <- cattle %>%
  group_by(week) %>%
  summarize(cattle_count = n())


###############
# combine date dataframes and export

HPAI_detections_week <- wildbird_count_week_wide %>%
  full_join(domes_count_week, by = "week") %>%
  full_join(mammals_count_week, by = "week") %>%
  full_join(cattle_count_week, by = "week") %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  mutate(galliformes = galliformes + domestic_count) %>%
  mutate(mammal = mammal_count + cattle_count) %>%
  select(-domestic_count) %>%
  rowwise() %>%
  ungroup()

write.csv(HPAI_detections_week, file = "HPAI_detections_week.csv", row.names = FALSE)




################
# visualize the number of seqs by order (sep anser into duck and goose)
# vissualize the number of detections by order
# rolling averages of =7 days and 3 weeks

colors <- c(
  accipitriformes   = "#0072B2",
  charadriiformes   = "#009E73",
  duck  = "#56B4E9",
  galliformes   = "#D55E00",
  goose = "#CC79A7",
  mammal = "#999933",
  passeriformes   = "#E69F00",
  strigiformes  = "#882255"
)




plot_df_week <- HPAI_detections_week %>%
  select(week, all_of(names(colors))) %>%
  pivot_longer(
    cols = -week,
    names_to = "group",
    values_to = "count") %>%
  arrange(group, week) %>%
  group_by(group) %>%
  mutate(
    rolling_avg = slide_dbl(count, mean, .before = 1, .after = 1, .complete = FALSE)) %>%
  ungroup()


ggplot(plot_df_week, aes(x = week, y = rolling_avg, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  labs(
    title = "Weekly Moving Average of HPAI Detections",
    x = "Week",
    y = "Moving Average (3-week)",
    color = "Group") +
  theme_minimal(base_size = 14)



########
# sequences over time

colors <- c(
  accipitriformes = "#0072B2",
  charadriiformes = "#009E73",
  duck            = "#56B4E9",
  galliformes     = "#D55E00",
  goose           = "#CC79A7",
  `nonhuman-mammal` = "#999933",
  passeriformes   = "#E69F00",
  strigiformes    = "#882255"
)


meta_daily <- meta %>%
  mutate(date = as.Date(date)) %>%
  filter(!is.na(date) & !is.na(order_dg)) %>% 
  count(date, order_dg, name = "count") %>%
  arrange(order_dg, date)

meta_daily <- meta_daily %>%
  group_by(order_dg) %>%
  mutate(
    rolling_avg = slide_dbl(count, mean, .before = 3, .after = 3, .complete = FALSE)) %>%
  ungroup()

ggplot(meta_daily, aes(x = date, y = rolling_avg, color = order_dg)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  labs(
    title = "Number of Sequences by host",
    x = "Date",
    y = "7-Day Rolling Average",
    color = "Group") +
  theme_minimal(base_size = 14)




meta_weekly <- meta %>%
  mutate(
    date = as.Date(date),
    week = floor_date(date, "week", week_start = 7)) %>%
  filter(!is.na(date) & !is.na(order_dg)) %>%
  count(week, order_dg, name = "count") %>%
  arrange(order_dg, week)


meta_weekly <- meta_weekly %>%
  group_by(order_dg) %>%
  mutate(
    rolling_avg = slide_dbl(count, mean, .before = 1, .after = 1, .complete = FALSE)) %>%
  ungroup()


ggplot(meta_weekly, aes(x = week, y = rolling_avg, color = order_dg)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  labs(
    title = "Number of Sequences by host",
    x = "Week",
    y = "3-Week Rolling Average",
    color = "Group") +
  theme_minimal(base_size = 14)


# combined plot
HPAI_week <- HPAI_detections_week %>%
  select(week, all_of(names(colors))) %>%
  pivot_longer(
    cols = -week,
    names_to = "group",
    values_to = "count") %>%
  arrange(group, week) %>%
  group_by(group) %>%
  mutate(
    rolling_avg = slide_dbl(count, mean, .before = 1, .after = 1, .complete = FALSE)) %>%
  ungroup()


meta_week <- meta %>%
  mutate(
    date = as.Date(date),
    week = floor_date(date, "week", week_start = 7)) %>%
  filter(!is.na(week) & !is.na(order_dg)) %>%
  count(week, order_dg, name = "count") %>%
  arrange(order_dg, week) %>%
  group_by(order_dg) %>%
  mutate(
    rolling_avg = slide_dbl(count, mean, .before = 1, .after = 1, .complete = FALSE)) %>%
  ungroup()


x_range <- range(c(meta_week$week, meta_week$week), na.rm = TRUE)

colors <- c(
  accipitriformes   = "#0072B2",
  charadriiformes   = "#009E73",
  duck  = "#56B4E9",
  galliformes   = "#D55E00",
  goose = "#CC79A7",
  mammal = "#999933",
  passeriformes   = "#E69F00",
  strigiformes  = "#882255"
)

p1 <- ggplot(HPAI_week, aes(x = week, y = rolling_avg, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  scale_x_date(limits = x_range, date_breaks = "4 month", date_labels = "%b\n%Y") +
  labs(
    title = "HPAI Detections (Weekly Rolling Avg)",
    x = NULL,
    y = "Count",
    color = "Group") +
  theme_minimal(base_size = 18)


colors <- c(
  accipitriformes = "#0072B2",
  charadriiformes = "#009E73",
  duck            = "#56B4E9",
  galliformes     = "#D55E00",
  goose           = "#CC79A7",
  `nonhuman-mammal` = "#999933",
  passeriformes   = "#E69F00",
  strigiformes    = "#882255"
)

p2 <- ggplot(meta_week, aes(x = week, y = rolling_avg, color = order_dg)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  scale_x_date(limits = x_range, date_breaks = "4 month", date_labels = "%b\n%Y") +
  labs(
    title = "Sequences by Host",
    x = "Date",
    y = "Count",
    color = "Group") +
  theme_minimal(base_size = 18)

combined_plot <- p1 / p2 + plot_layout(heights = c(1, 1))
combined_plot

ggsave("detections_sequences_HPAI_overtime_wcanada.pdf", plot = combined_plot, height = 8.5, width = 11, units = "in")



combined <- full_join(meta_monthly, wild_monthly, by = "month") %>%
  arrange(month) %>%
  replace_na(list(seq_count = 0, wild_count = 0)) %>%
  mutate(fraction = ifelse(wild_count == 0, NA, seq_count / wild_count))

combined <- combined %>%
  mutate(fraction_interp = na.approx(fraction, x = month, na.rm = FALSE))

plote <- ggplot(combined, aes(x = month, y = fraction_interp)) +
  geom_line(color = "darkgreen", size = 1.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Month",
    y = "Fraction (seq / wildbirds)",
    title = "Number of seqs / Number detections") +
  theme_minimal(base_size = 18)

plote

ggsave("seqs_vs_detections_overtimeplot.pdf", plot = plote, height = 8.5, width = 11, units = "in")


