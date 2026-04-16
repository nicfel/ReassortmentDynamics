#!/usr/bin/env Rscript

library(dplyr)

# Read the HPAI_LPAI.csv file
clade_file <- read.csv("./tables/HPAI_LPAI.csv", stringsAsFactors = FALSE, header = TRUE)

# Clean taxa names (remove quotes and parentheses)
clade_file$taxa <- gsub("'", "", clade_file$taxa)
clade_file$taxa <- gsub("\\(", "", clade_file$taxa)
clade_file$taxa <- gsub("\\)", "", clade_file$taxa)

# Extract location from taxa name (third element after splitting by "/")
clade_file$location <- sapply(clade_file$taxa, function(x) {
  parts <- strsplit(x, "/")[[1]]
  if(length(parts) >= 3) {
    return(parts[3])
  } else {
    return("Unknown")
  }
})

# Define flyway mappings based on US Fish and Wildlife Service flyway boundaries
# Pacific Flyway: Alaska, Hawaii, Washington, Oregon, California, Nevada, Idaho, Montana, Wyoming, Utah, Colorado, New Mexico, Arizona
# Central Flyway: North Dakota, South Dakota, Nebraska, Kansas, Oklahoma, Texas, Minnesota, Iowa, Missouri, Arkansas, Louisiana
# Mississippi Flyway: Montana (eastern part), North Dakota, Minnesota, Wisconsin, Michigan, Illinois, Indiana, Ohio, Kentucky, Tennessee, Mississippi, Alabama, Louisiana, Arkansas, Missouri, Iowa
# Atlantic Flyway: Maine, New Hampshire, Vermont, Massachusetts, Rhode Island, Connecticut, New York, New Jersey, Pennsylvania, Delaware, Maryland, Virginia, West Virginia, North Carolina, South Carolina, Georgia, Florida

# Note: Some states span multiple flyways, but we'll use the primary flyway

flyway_mapping <- list(
  # Pacific Flyway
  "Pacific" = c("Alaska", "Hawaii", "Washington", "Oregon", "California", "Nevada", "Idaho", 
                "Montana", "Wyoming", "Utah", "Colorado", "New_Mexico", "Arizona", "BC", "AB", "SK", "MB"),
  
  # Central Flyway  
  "Central" = c("North_Dakota", "South_Dakota", "Nebraska", "Kansas", "Oklahoma", "Texas", 
                "Minnesota", "Iowa", "Missouri", "Arkansas", "Louisiana", "New_Mexico"),
  
  # Mississippi Flyway
  "Mississippi" = c("Montana", "North_Dakota", "Minnesota", "Wisconsin", "Michigan", "Illinois", 
                    "Indiana", "Ohio", "Kentucky", "Tennessee", "Mississippi", "Alabama", 
                    "Louisiana", "Arkansas", "Missouri", "Iowa"),
  
  # Atlantic Flyway
  "Atlantic" = c("Maine", "New_Hampshire", "Vermont", "Massachusetts", "Rhode_Island", "Connecticut", 
                 "New_York", "New_Jersey", "Pennsylvania", "Delaware", "Maryland", "Virginia", 
                 "West_Virginia", "North_Carolina", "South_Carolina", "Georgia", "Florida", 
                 "NL", "NS", "NB", "PE", "QC", "ON", "NU"),
  
  # International locations
  "International" = c("Puebla", "Durango", "USA", "90")
)

# Function to assign flyway
assign_flyway <- function(location) {
  location <- gsub(" ", "_", location)
  
  for(flyway in names(flyway_mapping)) {
    if(location %in% flyway_mapping[[flyway]]) {
      return(flyway)
    }
  }
  return("Unknown")
}

# Assign flyways
clade_file$flyway <- sapply(clade_file$location, assign_flyway)

# Handle special cases and overlaps
# Some states span multiple flyways - we'll use more specific assignments
special_cases <- list(
  "Minnesota" = "Central",  # Primarily Central Flyway
  "Iowa" = "Central",       # Primarily Central Flyway  
  "Missouri" = "Central",   # Primarily Central Flyway
  "Arkansas" = "Central",   # Primarily Central Flyway
  "Louisiana" = "Central",  # Primarily Central Flyway
  "Montana" = "Pacific",    # Primarily Pacific Flyway
  "North_Dakota" = "Central", # Primarily Central Flyway
  "New_Mexico" = "Central"  # Primarily Central Flyway
)

# Apply special cases
for(loc in names(special_cases)) {
  clade_file$flyway[clade_file$location == loc] <- special_cases[[loc]]
}

# Summary statistics
cat("=== FLYWAY ASSIGNMENT SUMMARY ===\n")
flyway_summary <- table(clade_file$flyway)
print(flyway_summary)

cat("\n=== LOCATION SUMMARY BY FLYWAY ===\n")
location_summary <- clade_file %>%
  group_by(flyway, location) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(flyway, desc(count))

print(location_summary)

# Check for any unassigned locations
unassigned <- clade_file[clade_file$flyway == "Unknown", ]
if(nrow(unassigned) > 0) {
  cat("\n=== UNASSIGNED LOCATIONS ===\n")
  print(unique(unassigned$location))
}

# Create a summary by status and flyway
cat("\n=== STATUS BY FLYWAY ===\n")
status_flyway_summary <- table(clade_file$status, clade_file$flyway)
print(status_flyway_summary)

# Save the enhanced dataset
write.csv(clade_file, "./tables/HPAI_LPAI_with_flyways.csv", row.names = FALSE)
cat("\nEnhanced dataset saved to: ./tables/HPAI_LPAI_with_flyways.csv\n")

# Create a flyway mapping reference file
flyway_reference <- data.frame(
  location = unlist(lapply(names(flyway_mapping), function(flyway) {
    rep(flyway, length(flyway_mapping[[flyway]]))
  })),
  state_province = unlist(flyway_mapping)
)

write.csv(flyway_reference, "./tables/flyway_mapping_reference.csv", row.names = FALSE)
cat("Flyway mapping reference saved to: ./tables/flyway_mapping_reference.csv\n")
