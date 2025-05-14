# Load required libraries
library(dplyr)
library(readr)

# ---- USER INPUTS ----
input_file <- "your_data.csv"  # Path to your CSV file
target_individuals <- 1652     # Number of individuals needed per replicate
max_volume <- 9                # Max allowable total volume in mL
supplement_density <- 630      # Density of the secondary source (individuals/mL)

# ---- LOAD DATA ----
data <- read_csv(input_file)

# ---- PROCESSING ----
output <- data %>%
  mutate(
    required_volume = target_individuals / count_ml,
    needs_supplement = required_volume > max_volume,
    volume_primary = ifelse(
      needs_supplement,
      max_volume - (target_individuals - max_volume * count_ml) / (supplement_density - count_ml),
      required_volume
    ),
    volume_secondary = ifelse(
      needs_supplement,
      max_volume - volume_primary,
      0
    ),
    total_individuals = volume_primary * count_ml + volume_secondary * supplement_density,
    total_volume = volume_primary + volume_secondary
  )

# ---- OUTPUT ----
# Print or save the result
print(output)
# write_csv(output, "supplemented_cultures.csv")  # Uncomment to save output
