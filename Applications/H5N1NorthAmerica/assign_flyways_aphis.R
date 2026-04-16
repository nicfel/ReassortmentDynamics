#!/usr/bin/env Rscript

library(dplyr)
library(lubridate)
library(ggplot2)
library(maps)
library(viridis)

# Read the APHIS Wild Bird Avian Influenza Surveillance Dashboard data
cat("Reading APHIS data...\n")
aphis_data <- read.csv("./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv", stringsAsFactors = FALSE, header = TRUE)

cat("Dataset dimensions:", dim(aphis_data), "\n")
cat("Column names:", colnames(aphis_data), "\n")

# Clean the State column (remove any leading/trailing whitespace)
aphis_data$State <- trimws(aphis_data$State)

# Get unique states
unique_states <- sort(unique(aphis_data$State))
cat("\nUnique states in dataset:", length(unique_states), "\n")
print(unique_states)

# Define comprehensive flyway mappings based on US Fish and Wildlife Service flyway boundaries
# Note: Some states span multiple flyways, but we'll use the primary flyway for each state

flyway_mapping <- list(
  # Pacific Flyway
  "Pacific" = c("Alaska", "Hawaii", "Washington", "Oregon", "California", "Nevada", "Idaho", 
                "Arizona", "Utah"),
  
  # Central Flyway  
  "Central" = c("North Dakota", "South Dakota", "Nebraska", "Kansas", "Oklahoma", "Texas", 
                "Montana", "Wyoming", "Colorado", "New Mexico"),
  
  # Mississippi Flyway
  "Mississippi" = c("Minnesota", "Wisconsin", "Michigan", "Illinois", 
                    "Indiana", "Ohio", "Kentucky", "Tennessee", "Mississippi", "Alabama", 
                    "Louisiana", "Arkansas", "Missouri", "Iowa"),
  
  # Atlantic Flyway
  "Atlantic" = c("Maine", "New Hampshire", "Vermont", "Massachusetts", "Rhode Island", "Connecticut", 
                 "New York", "New Jersey", "Pennsylvania", "Delaware", "Maryland", "Virginia", 
                 "West Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", 
                 "Puerto Rico", "Virgin Islands")
)

# Function to assign flyway
assign_flyway <- function(state) {
  state <- trimws(state)
  
  for(flyway in names(flyway_mapping)) {
    if(state %in% flyway_mapping[[flyway]]) {
      return(flyway)
    }
  }
  return("Unknown")
}

# Assign flyways to all records
cat("\nAssigning flyways...\n")
aphis_data$flyway <- sapply(aphis_data$State, assign_flyway)

# Handle special cases where states span multiple flyways using county-level information
# Define county-level flyway assignments for states that span multiple flyways

# Based on official US Fish and Wildlife Service flyway boundaries
# Multi-flyway states with county-level assignments
county_flyway_mapping <- list(
  # 1. Texas: Central (western/central) & Mississippi (eastern third)
  "Texas" = list(
    "Central" = c("Andrews", "Armstrong", "Bailey", "Baylor", "Borden", "Brewster", "Briscoe", 
                  "Brown", "Callahan", "Carson", "Castro", "Childress", "Cochran", "Coke", 
                  "Coleman", "Collingsworth", "Concho", "Crane", "Crosby", "Culberson", "Dallam", 
                  "Dawson", "Deaf Smith", "Dickens", "Donley", "Eastland", "Ector", "El Paso", 
                  "Fisher", "Floyd", "Foard", "Gaines", "Garza", "Glasscock", "Gray", "Hale", 
                  "Hansford", "Hardeman", "Hartley", "Haskell", "Hemphill", "Hockley", "Howard", 
                  "Hudspeth", "Hutchinson", "Irion", "Jeff Davis", "Jones", "Kent", "King", 
                  "Knox", "Lamb", "Loving", "Lubbock", "Lynn", "Martin", "Mason", "McCulloch", 
                  "Midland", "Mitchell", "Moore", "Motley", "Nolan", "Pecos", "Potter", "Presidio", 
                  "Reagan", "Reeves", "Runnels", "Schleicher", "Scurry", "Shackelford", "Sherman", 
                  "Sterling", "Stonewall", "Sutton", "Swisher", "Taylor", "Terrell", "Terry", 
                  "Throckmorton", "Tom Green", "Upton", "Ward", "Wheeler", "Winkler", "Yoakum"),
    "Mississippi" = c("Angelina", "Atascosa", "Austin", "Bandera", "Bastrop", "Bee", "Bell", 
                      "Bexar", "Blanco", "Bosque", "Bowie", "Brazoria", "Brazos", "Burleson", 
                      "Burnet", "Caldwell", "Calhoun", "Cameron", "Camp", "Cass", "Chambers", 
                      "Cherokee", "Clay", "Colorado", "Comal", "Comanche", "Cooke", "Coryell", 
                      "Dallas", "Delta", "Denton", "DeWitt", "Duval", "Ellis", "Erath", "Falls", 
                      "Fannin", "Fayette", "Fort Bend", "Franklin", "Freestone", "Galveston", 
                      "Gillespie", "Goliad", "Gonzales", "Grayson", "Gregg", "Grimes", "Guadalupe", 
                      "Hardin", "Harris", "Harrison", "Hays", "Henderson", "Hidalgo", "Hill", 
                      "Hood", "Hopkins", "Houston", "Hunt", "Jackson", "Jasper", "Jefferson", 
                      "Jim Hogg", "Jim Wells", "Johnson", "Karnes", "Kaufman", "Kendall", "Kerr", 
                      "Kimble", "Kleberg", "Lamar", "Lampasas", "Lavaca", "Lee", "Leon", "Liberty", 
                      "Limestone", "Live Oak", "Llano", "Madison", "Marion", "Matagorda", "Maverick", 
                      "McLennan", "McMullen", "Medina", "Menard", "Milam", "Montgomery", "Morris", 
                      "Nacogdoches", "Navarro", "Newton", "Orange", "Palo Pinto", "Panola", "Parker", 
                      "Polk", "Rains", "Red River", "Refugio", "Robertson", "Rockwall", "Rusk", 
                      "Sabine", "San Augustine", "San Jacinto", "San Patricio", "San Saba", 
                      "Shelby", "Smith", "Somervell", "Starr", "Tarrant", "Titus", "Travis", 
                      "Trinity", "Tyler", "Upshur", "Uvalde", "Val Verde", "Van Zandt", "Victoria", 
                      "Walker", "Waller", "Washington", "Webb", "Wharton", "Wichita", "Wilbarger", 
                      "Willacy", "Williamson", "Wilson", "Wise", "Wood", "Young", "Zapata", "Zavala")
  ),
  
  # 2. Oklahoma: Central (western/central) & Mississippi (eastern)
  "Oklahoma" = list(
    "Central" = c("Alfalfa", "Beaver", "Beckham", "Blaine", "Caddo", "Canadian", "Cimarron", 
                  "Cleveland", "Comanche", "Custer", "Dewey", "Ellis", "Garfield", "Grant", 
                  "Greer", "Harmon", "Harper", "Jackson", "Jefferson", "Kay", "Kingfisher", 
                  "Kiowa", "Logan", "Major", "Noble", "Oklahoma", "Payne", "Pottawatomie", 
                  "Roger Mills", "Seminole", "Texas", "Tillman", "Washita", "Woods", "Woodward"),
    "Mississippi" = c("Adair", "Atoka", "Bryan", "Cherokee", "Choctaw", "Coal", "Craig", 
                      "Creek", "Delaware", "Garvin", "Grady", "Haskell", "Hughes", "Johnston", 
                      "Latimer", "Le Flore", "Lincoln", "Love", "McClain", "McCurtain", "McIntosh", 
                      "Marshall", "Mayes", "Muskogee", "Nowata", "Okfuskee", "Okmulgee", "Osage", 
                      "Ottawa", "Pittsburg", "Pontotoc", "Pushmataha", "Rogers", "Sequoyah", 
                      "Stephens", "Tulsa", "Wagoner", "Washington")
  ),
  
  # 3. Kansas: Central (western) & Mississippi (eastern)
  "Kansas" = list(
    "Central" = c("Barber", "Barton", "Cheyenne", "Clark", "Comanche", "Decatur", "Edwards", 
                  "Ellis", "Finney", "Ford", "Gove", "Graham", "Grant", "Gray", "Greeley", 
                  "Hamilton", "Harper", "Haskell", "Hodgeman", "Jewell", "Kearny", "Kingman", 
                  "Kiowa", "Lane", "Logan", "Meade", "Mitchell", "Morton", "Ness", "Norton", 
                  "Pawnee", "Phillips", "Rawlins", "Reno", "Rooks", "Rush", "Russell", "Scott", 
                  "Seward", "Sheridan", "Sherman", "Stafford", "Stanton", "Stevens", "Thomas", 
                  "Trego", "Wallace", "Wichita", "Wyandotte"),
    "Mississippi" = c("Allen", "Anderson", "Atchison", "Bourbon", "Brown", "Butler", "Chase", 
                      "Chautauqua", "Clay", "Cloud", "Coffey", "Cowley", "Crawford", "Dickinson", 
                      "Doniphan", "Douglas", "Elk", "Franklin", "Geary", "Jackson", "Jefferson", 
                      "Johnson", "Labette", "Leavenworth", "Linn", "Lyon", "Marion", "Marshall", 
                      "McPherson", "Miami", "Montgomery", "Morris", "Nemaha", "Neosho", "Osage", 
                      "Ottawa", "Pottawatomie", "Pratt", "Republic", "Riley", "Saline", "Sedgwick", 
                      "Shawnee", "Sumner", "Washington", "Wilson", "Woodson")
  ),
  
  # 4. Nebraska: Central (western) & Mississippi (eastern)
  "Nebraska" = list(
    "Central" = c("Adams", "Arthur", "Banner", "Blaine", "Box Butte", "Boyd", "Brown", "Buffalo", 
                  "Butler", "Cass", "Cherry", "Cheyenne", "Cuming", "Custer", "Dakota", "Dawes", 
                  "Dawson", "Deuel", "Dixon", "Dundy", "Fillmore", "Frontier", "Furnas", "Garden", 
                  "Garfield", "Gosper", "Grant", "Greeley", "Hall", "Hamilton", "Harlan", "Hayes", 
                  "Hitchcock", "Holt", "Hooker", "Howard", "Jefferson", "Johnson", "Kearney", 
                  "Keith", "Keya Paha", "Kimball", "Knox", "Lancaster", "Lincoln", "Logan", 
                  "Loup", "McPherson", "Madison", "Merrick", "Morrill", "Nance", "Nuckolls", 
                  "Otoe", "Pawnee", "Perkins", "Phelps", "Pierce", "Platte", "Polk", "Red Willow", 
                  "Richardson", "Rock", "Saline", "Sarpy", "Scotts Bluff", "Seward", "Sheridan", 
                  "Sherman", "Sioux", "Stanton", "Thayer", "Thomas", "Thurston", "Valley", 
                  "Washington", "Wayne", "Webster", "Wheeler", "York"),
    "Mississippi" = c("Antelope", "Boone", "Burt", "Cedar", "Clay", "Colfax", "Dodge", "Douglas", 
                      "Gage", "Greeley", "Madison", "Nance", "Nemaha", "Otoe", "Platte", "Polk", 
                      "Richardson", "Saline", "Sarpy", "Saunders", "Seward", "Thurston", "Washington", 
                      "Wayne", "York")
  ),
  
  # 5. North Dakota: Central (western) & Mississippi (eastern)
  "North Dakota" = list(
    "Central" = c("Adams", "Billings", "Bowman", "Burke", "Divide", "Dunn", "Golden Valley", 
                  "Grant", "Hettinger", "McKenzie", "Mercer", "Morton", "Mountrail", "Oliver", 
                  "Renville", "Slope", "Stark", "Ward", "Williams", "Bottineau", "McHenry", 
                  "Pierce", "Rolette", "Towner", "Burleigh", "Emmons", "Kidder", "McLean", 
                  "Sheridan", "Sioux", "Wells"),
    "Mississippi" = c("Cass", "Grand Forks", "Nelson", "Pembina", "Ramsey", "Richland", 
                      "Sargent", "Steele", "Traill", "Walsh")
  ),
  
  # 6. South Dakota: Central (western) & Mississippi (eastern)
  "South Dakota" = list(
    "Central" = c("Beadle", "Bennett", "Bon Homme", "Brookings", "Brown", "Brule", "Buffalo", 
                  "Campbell", "Charles Mix", "Clark", "Clay", "Codington", "Corson", "Custer", 
                  "Davison", "Day", "Deuel", "Dewey", "Douglas", "Edmunds", "Fall River", 
                  "Faulk", "Grant", "Gregory", "Haakon", "Hamlin", "Hand", "Hanson", "Harding", 
                  "Hughes", "Hutchinson", "Hyde", "Jackson", "Jerauld", "Jones", "Kingsbury", 
                  "Lake", "Lawrence", "Lincoln", "Lyman", "Marshall", "McCook", "McPherson", 
                  "Meade", "Mellette", "Miner", "Minnehaha", "Moody", "Oglala Lakota", "Pennington", 
                  "Perkins", "Potter", "Roberts", "Sanborn", "Shannon", "Spink", "Stanley", 
                  "Sully", "Todd", "Tripp", "Turner", "Union", "Walworth", "Yankton", "Ziebach"),
    "Mississippi" = c("Brookings", "Clay", "Codington", "Hamlin", "Kingsbury", "Lake", "Lincoln", 
                      "Marshall", "McCook", "Minnehaha", "Moody", "Roberts", "Turner", "Union", "Yankton")
  ),
  
  
  # 11. Louisiana: Entirely Mississippi Flyway (all parishes)
  "Louisiana" = list(
    "Mississippi" = c("Acadia", "Allen", "Ascension", "Assumption", "Avoyelles", "Beauregard", 
                      "Bienville", "Bossier", "Caddo", "Calcasieu", "Caldwell", "Cameron", 
                      "Catahoula", "Claiborne", "Concordia", "De Soto", "East Baton Rouge", 
                      "East Carroll", "Evangeline", "Franklin", "Grant", "Iberia", "Iberville", 
                      "Jackson", "Jefferson", "Jefferson Davis", "LaSalle", "Lafayette", 
                      "Lafourche", "Lincoln", "Livingston", "Madison", "Morehouse", "Natchitoches", 
                      "Orleans", "Ouachita", "Plaquemines", "Pointe Coupee", "Rapides", 
                      "Red River", "Richland", "Sabine", "St. Bernard", "St. Charles", 
                      "St. Helena", "St. James", "St. John the Baptist", "St. Landry", 
                      "St. Martin", "St. Mary", "St. Tammany", "Tangipahoa", "Tensas", 
                      "Terrebonne", "Union", "Vermilion", "Vernon", "Washington", "Webster", 
                      "West Baton Rouge", "West Carroll", "West Feliciana", "Winn")
  ),
  
  
  # 16. Alaska: Pacific (western/southern coast) & Central (interior)
  "Alaska" = list(
    "Pacific" = c("Aleutians East", "Aleutians West", "Anchorage", "Bethel", "Bristol Bay", 
                  "Dillingham", "Haines", "Juneau", "Kenai Peninsula", "Ketchikan Gateway", 
                  "Kodiak Island", "Lake and Peninsula", "Matanuska-Susitna", "Nome", 
                  "North Slope", "Northwest Arctic", "Prince of Wales-Hyder", "Sitka", 
                  "Skagway", "Southeast Fairbanks", "Valdez-Cordova", "Wade Hampton", "Wrangell", 
                  "Yakutat", "Yukon-Koyukuk"),
    "Central" = c("Fairbanks North Star", "Southeast Fairbanks", "Yukon-Koyukuk")
  )
)

# Function to assign flyway based on county
assign_county_flyway <- function(state, county) {
  state <- trimws(state)
  county <- trimws(county)
  
  # Check if this state has county-level mapping
  if(state %in% names(county_flyway_mapping)) {
    state_mapping <- county_flyway_mapping[[state]]
    
    # Check each flyway for this county
    for(flyway in names(state_mapping)) {
      if(county %in% state_mapping[[flyway]]) {
        return(flyway)
      }
    }
  }
  
  # If no county mapping found, return the default state flyway
  return(assign_flyway(state))
}

# Apply county-level flyway assignments for special cases
cat("\nApplying county-level flyway assignments...\n")
special_case_states <- names(county_flyway_mapping)

for(state in special_case_states) {
  state_mask <- aphis_data$State == state
  if(sum(state_mask) > 0) {
    aphis_data$flyway[state_mask] <- mapply(assign_county_flyway, 
                                           aphis_data$State[state_mask], 
                                           aphis_data$County[state_mask])
  }
}

# Summary statistics
cat("\n=== FLYWAY ASSIGNMENT SUMMARY ===\n")
flyway_summary <- table(aphis_data$flyway)
print(flyway_summary)

cat("\n=== STATE SUMMARY BY FLYWAY ===\n")
state_flyway_summary <- aphis_data %>%
  group_by(flyway, State) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(flyway, desc(count))

print(state_flyway_summary)

# Show county-level assignments for special case states
cat("\n=== COUNTY-LEVEL FLYWAY ASSIGNMENTS ===\n")
for(state in special_case_states) {
  state_data <- aphis_data[aphis_data$State == state, ]
  if(nrow(state_data) > 0) {
    county_summary <- state_data %>%
      group_by(County, flyway) %>%
      summarise(count = n(), .groups = 'drop') %>%
      arrange(County, flyway)
    
    cat("\n", state, ":\n")
    print(county_summary)
  }
}

# Check for any unassigned states
unassigned <- aphis_data[aphis_data$flyway == "Unknown", ]
if(nrow(unassigned) > 0) {
  cat("\n=== UNASSIGNED STATES ===\n")
  unassigned_states <- unique(unassigned$State)
  print(unassigned_states)
}

# Create summary by detection status and flyway
cat("\n=== DETECTION STATUS BY FLYWAY ===\n")
# Convert date to proper format
aphis_data$Date_Collected <- as.Date(aphis_data$Date_Collected, format="%Y-%m-%d")

# Summary by Final_H5 detection
h5_summary <- table(aphis_data$Final_H5, aphis_data$flyway, useNA = "ifany")
print(h5_summary)

# Summary by Final_IAV detection
iav_summary <- table(aphis_data$Final_IAV, aphis_data$flyway, useNA = "ifany")
print(iav_summary)

# Summary by pathogenicity
pathogenicity_summary <- table(aphis_data$Final_Pathogenicity, aphis_data$flyway, useNA = "ifany")
print(pathogenicity_summary)

# Create yearly summary by flyway
cat("\n=== YEARLY SUMMARY BY FLYWAY ===\n")
yearly_flyway <- aphis_data %>%
  group_by(year, flyway) %>%
  summarise(
    total_samples = n(),
    h5_detected = sum(Final_H5 == "Detected", na.rm = TRUE),
    iav_detected = sum(Final_IAV == "Detected", na.rm = TRUE),
    hpai_detected = sum(Final_Pathogenicity == "High Path AI", na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(year, flyway)

print(yearly_flyway)

# Save the enhanced dataset
write.csv(aphis_data, "./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard_with_flyways.csv", row.names = FALSE)
cat("\nEnhanced dataset saved to: ./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard_with_flyways.csv\n")

# Create a flyway mapping reference file
flyway_reference <- data.frame(
  flyway = unlist(lapply(names(flyway_mapping), function(flyway) {
    rep(flyway, length(flyway_mapping[[flyway]]))
  })),
  state = unlist(flyway_mapping),
  stringsAsFactors = FALSE
)

write.csv(flyway_reference, "./tables/flyway_mapping_reference_aphis.csv", row.names = FALSE)
cat("Flyway mapping reference saved to: ./tables/flyway_mapping_reference_aphis.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Total records processed:", nrow(aphis_data), "\n")
cat("Date range:", min(aphis_data$Date_Collected, na.rm = TRUE), "to", max(aphis_data$Date_Collected, na.rm = TRUE), "\n")
cat("States covered:", length(unique(aphis_data$State)), "\n")
cat("Flyways assigned:", length(unique(aphis_data$flyway)), "\n")

# Create a US map visualization of flyway assignments
cat("\n=== CREATING FLYWAY MAP VISUALIZATION ===\n")

# Create a data frame with state-level flyway assignments
state_flyway_summary <- aphis_data %>%
  group_by(State, flyway) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(State) %>%
  slice_max(count, n = 1) %>%  # Get the dominant flyway for each state
  ungroup()

# Get US states map data
us_states <- map_data("state")
us_states$region <- tools::toTitleCase(us_states$region)

# Create a mapping of state names to flyways
# Handle special cases where states might be split
state_to_flyway <- setNames(state_flyway_summary$flyway, state_flyway_summary$State)

# Add flyway information to the map data
us_states$flyway <- state_to_flyway[us_states$region]
us_states$flyway[is.na(us_states$flyway)] <- "Unknown"

# Define colors for each flyway
flyway_colors <- c(
  "Pacific" = "#E31A1C",      # Red
  "Central" = "#1F78B4",      # Blue  
  "Mississippi" = "#33A02C",  # Green
  "Atlantic" = "#FF7F00",     # Orange
  "Unknown" = "#B15928"       # Brown
)

# Create the map
flyway_map <- ggplot(us_states, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = flyway), color = "white", size = 0.2) +
  scale_fill_manual(
    values = flyway_colors,
    name = "Flyway",
    labels = c("Atlantic", "Central", "Mississippi", "Pacific", "Unknown")
  ) +
  coord_fixed(1.3) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    title = "US Flyway Assignments",
    subtitle = "Based on APHIS Wild Bird Surveillance Data"
  )

# Save the map
ggsave("./tables/US_Flyway_Map.png", flyway_map, width = 12, height = 8, dpi = 300)
cat("Flyway map saved to: ./tables/US_Flyway_Map.png\n")

# Create a county-level map for multi-flyway states (if possible)
cat("\nCreating detailed county-level map for multi-flyway states...\n")

# Get county map data
us_counties <- map_data("county")

# Create a comprehensive county-to-flyway mapping
county_to_flyway <- data.frame()

# Process each multi-flyway state
for(state in names(county_flyway_mapping)) {
  state_mapping <- county_flyway_mapping[[state]]
  state_lower <- tolower(state)
  
  for(flyway in names(state_mapping)) {
    counties <- state_mapping[[flyway]]
    counties_lower <- tolower(counties)
    
    # Add to mapping
    county_df <- data.frame(
      state = rep(state_lower, length(counties_lower)),
      county = counties_lower,
      flyway = rep(flyway, length(counties_lower)),
      stringsAsFactors = FALSE
    )
    county_to_flyway <- rbind(county_to_flyway, county_df)
  }
}

# Merge county data with flyway assignments
us_counties$flyway <- NA
for(i in 1:nrow(county_to_flyway)) {
  mask <- us_counties$region == county_to_flyway$state[i] & 
          us_counties$subregion == county_to_flyway$county[i]
  us_counties$flyway[mask] <- county_to_flyway$flyway[i]
}

# Fill in single-flyway states
single_flyway_states <- state_flyway_summary[!state_flyway_summary$State %in% names(county_flyway_mapping), ]
for(i in 1:nrow(single_flyway_states)) {
  state_lower <- tolower(single_flyway_states$State[i])
  mask <- us_counties$region == state_lower & is.na(us_counties$flyway)
  us_counties$flyway[mask] <- single_flyway_states$flyway[i]
}

# Special cases: These states should be entirely assigned to their respective flyways
# If any counties are still unknown, assign them appropriately
mississippi_states <- c("louisiana", "minnesota", "iowa", "missouri", "arkansas")
for(state in mississippi_states) {
  state_mask <- us_counties$region == state & is.na(us_counties$flyway)
  us_counties$flyway[state_mask] <- "Mississippi"
}

central_states <- c("montana", "wyoming", "colorado", "new mexico")
for(state in central_states) {
  state_mask <- us_counties$region == state & is.na(us_counties$flyway)
  us_counties$flyway[state_mask] <- "Central"
}

# Function to assign unknown counties based on surrounding counties
assign_unknown_by_neighbors <- function(us_counties) {
  # Get all unknown counties
  unknown_mask <- is.na(us_counties$flyway)
  unknown_counties <- us_counties[unknown_mask, ]
  
  if(nrow(unknown_counties) == 0) return(us_counties)
  
  cat("Assigning", nrow(unknown_counties), "unknown counties based on neighbors...\n")
  
  for(i in 1:nrow(unknown_counties)) {
    county_state <- unknown_counties$region[i]
    
    # Get all counties in the same state
    state_counties <- us_counties[us_counties$region == county_state & !is.na(us_counties$flyway), ]
    
    if(nrow(state_counties) > 0) {
      # Use the most common flyway in the state
      flyway_counts <- table(state_counties$flyway)
      most_common_flyway <- names(flyway_counts)[which.max(flyway_counts)]
      
      # Assign to the most common flyway in the state
      county_idx <- which(us_counties$region == county_state & 
                         us_counties$subregion == unknown_counties$subregion[i])
      us_counties$flyway[county_idx] <- most_common_flyway
      
      cat("Assigned", unknown_counties$subregion[i], ",", county_state, "to", most_common_flyway, "flyway\n")
    }
  }
  
  return(us_counties)
}

# Assign unknown counties based on neighbors
us_counties <- assign_unknown_by_neighbors(us_counties)

# Set any remaining unknown counties (should be very few now)
us_counties$flyway[is.na(us_counties$flyway)] <- "Unknown"

# Create county-level map
county_flyway_map <- ggplot(us_counties, aes(x = long, y = lat, group = group)) +
  geom_polygon(aes(fill = flyway), color = "white", size = 0.1) +
  scale_fill_manual(
    values = flyway_colors,
    name = "Flyway",
    labels = c("Atlantic", "Central", "Mississippi", "Pacific", "Unknown")
  ) +
  coord_fixed(1.3) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    title = "US Flyway Assignments (County Level)",
    subtitle = "Detailed view showing county-level flyway boundaries"
  )

# Save the county-level map
ggsave("./tables/US_Flyway_Map_County_Level.png", county_flyway_map, width = 12, height = 8, dpi = 300)
cat("County-level flyway map saved to: ./tables/US_Flyway_Map_County_Level.png\n")

# Print summary of county-level assignments
cat("\n=== COUNTY-LEVEL FLYWAY ASSIGNMENT SUMMARY ===\n")
county_flyway_summary <- us_counties %>%
  group_by(flyway) %>%
  summarise(counties = n(), .groups = 'drop') %>%
  arrange(desc(counties))
print(county_flyway_summary)

# Debug: Show some counties that might have mapping issues
cat("\n=== DEBUGGING: COUNTIES WITH POTENTIAL MAPPING ISSUES ===\n")
# Check for counties that are still "Unknown"
unknown_counties <- us_counties[us_counties$flyway == "Unknown", ]
if(nrow(unknown_counties) > 0) {
  cat("Counties marked as Unknown (first 20):\n")
  unknown_sample <- unknown_counties %>%
    group_by(region, subregion) %>%
    summarise(count = n(), .groups = 'drop') %>%
    head(20)
  print(unknown_sample)
}

# Check Louisiana specifically since we know it should be all Mississippi
cat("\n=== LOUISIANA COUNTY CHECK ===\n")
la_counties <- us_counties[us_counties$region == "louisiana", ]
la_summary <- la_counties %>%
  group_by(flyway) %>%
  summarise(count = n(), .groups = 'drop')
print(la_summary)
