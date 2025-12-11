library(dplyr)
library(lubridate)
library(scales)
library(viridis)
library(zoo)
library(patchwork)
library(ggplot2)
library(tidyr)
library(readr)

setwd("~/Documents/HPAI/H5-reassortment/prev-reassortment/H5N1_detections/")

#### format all other AIV data from USDA dashboard.
#https://www.aphis.usda.gov/livestock-poultry-disease/avian/avian-influenza/wild-bird-surveillance-dashboard

AIVdash <- read.csv("APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv")
# fix date column
AIVdash <- AIVdash %>%
  mutate(Date_Clean = as.Date(sub(" .*", "", Date_Collected), format = "%m/%d/%y"))

# Replace all blank strings in the dataframe with NA
AIVdash[AIVdash == ""] <- NA
AIVdash$Final_IAV <- as.character(AIVdash$Final_IAV)

overall_detected_pct <- mean(AIVdash$Final_IAV == "Detected", na.rm = TRUE) * 100
cat("Overall Detected Percentage:", round(overall_detected_pct, 2), "%\n")

AIVdash <- AIVdash %>%
  mutate(weekly = floor_date(Date_Clean, unit = "week"))



AIVdash <- AIVdash %>%
  mutate(Final_Pathogenicity_Clean = ifelse(is.na(Final_Pathogenicity), "LPAI", Final_Pathogenicity))
#Create time bins
AIVdash <- AIVdash %>%
  mutate(week = floor_date(Date_Clean, unit = "week"))



########
# flyway based propotion postivity
flyways <- read_csv("flyway_regions.csv")
AIVdash <- AIVdash %>%
  left_join(flyways, by = c("State" = "location"))


# Summarise by flyway & week
flyway_weekly_ratio <- AIVdash %>%
  filter(!is.na(flyway)) %>%
  group_by(flyway, week) %>%
  summarise(
    high_path_count = sum(Final_Pathogenicity == "High Path AI" & Final_H5 == "Detected", na.rm = TRUE),
    h5_detected_count = sum(Final_H5 == "Detected", na.rm = TRUE),
    ratio = ifelse(h5_detected_count > 0, high_path_count / h5_detected_count, 0)) %>%
  ungroup()


##############
# test posit ivity plots by flyway
# Summarize for each flyway & week
flyway_weekly_pos <- AIVdash %>%
  filter(!is.na(flyway)) %>%
  group_by(flyway, week) %>%
  summarise(total_samples = n(),
    h5_detected = sum(
      Final_H5 == "Detected" & Final_Pathogenicity_Clean != "High Path AI",
      na.rm = TRUE),
    hpai_detected = sum(
      Final_Pathogenicity == "High Path AI",
      na.rm = TRUE),
    h5_pos_rate = ifelse(total_samples > 0, h5_detected / total_samples, 0),
    hpai_pos_rate = ifelse(total_samples > 0, hpai_detected / total_samples, 0)) %>%
  ungroup()

######## detections by flyway
flyway_weekly_counts <- AIVdash %>%
  filter(!is.na(flyway)) %>%
  mutate(week = floor_date(Date_Clean, unit = "week")) %>%
  group_by(flyway, week) %>%
  summarise(
    low_path_h5 = sum(Final_H5 == "Detected" & Final_Pathogenicity_Clean == "LPAI", na.rm = TRUE),
    high_path_h5 = sum(Final_Pathogenicity == "High Path AI", na.rm = TRUE)) %>%
  ungroup()

flyway_weekly_counts_long <- flyway_weekly_counts %>%
  pivot_longer(
    cols = c(low_path_h5, high_path_h5),
    names_to = "path_type",
    values_to = "count")

path_colors <- c(
  "low_path_h5" = "#E69F00",
  "high_path_h5" = "#56B4E9")

ggplot(flyway_weekly_counts_long, aes(x = week, y = count, fill = path_type)) +
  geom_area(position = "stack") +
  facet_wrap(~ flyway) +
  scale_fill_manual(
    values = path_colors,
    name = "Path Type",
    labels = c("Low Path H5", "High Path AI")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
  labs(
    title = "Weekly Stacked H5 Detections (Low Path vs. High Path) by Flyway",
    x = "Week",
    y = "Number of Detections") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(flyway_weekly_counts_long, aes(x = week, y = count, color = path_type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ flyway) +
  scale_color_manual(
    values = path_colors,
    name = "Path Type",
    labels = c("Low Path H5", "High Path AI")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") +
  labs(
    title = "Weekly Stacked H5 Detections (Low Path vs. High Path) by Flyway",
    x = "Week",
    y = "Number of Detections") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####### total detections by flyway LPAI vs HPAI

ggplot(flyway_weekly_counts_long, aes(x = week, y = count, fill = path_type)) +
  geom_area(position = "stack") +
  scale_fill_manual(
    values = path_colors,
    name = "Path Type",
    labels = c("Low Path H5", "High Path AI")) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  labs(
    title = "Weekly Stacked H5 Detections (Low Path vs. High Path) by Flyway",
    x = "Week",
    y = "Number of Detections") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

