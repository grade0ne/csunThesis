library(tidyverse)

Protists <- c("Bodo", "Colpoda", "Chrysomonads", "Colpidium", "Cryptomonads", "Cyclidium")

data <- read.csv("Data/PP_succession_censusdata_cleaned.csv")

data$Bodo <- as.integer(data$Bodo)

data <- data %>%
  filter(if_any(all_of(Protists), ~ . > 0)) %>%
  mutate(
        plantID = as.factor(paste(Plant.ID.Alpha, Plant.ID.Numeric)),
        cooccurr = Rotifers > 0
  ) %>%
  select(plantID, all_of(Protists), Rotifers, cooccurr)

# want abundance of each individual protist group with and without rotifer co-occurrance

graphdata <- data %>%
  pivot_longer(
    cols = all_of(Protists),
    names_to = "Protist"
  ) %>%
  mutate(
    occurrence = if_else(cooccurr, "together", "solo")
  ) %>%
  group_by(Protist, occurrence) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value) / sqrt(length(value)),
    Rotifers = Rotifers
  )

ggplot(graphdata, aes(x = Protist, y = mean, group = occurrence, fill = occurrence)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(stat = "identity", aes(ymax = mean + se, ymin = mean - se),
                position = position_dodge())

# want multiple regression with count as response and occurence and protist as predictors, plant ID as random factor

regData <- data %>%
  pivot_longer(
    cols = all_of(Protists),
    names_to = "Protist"
  ) %>%
  mutate(occurrence = if_else(cooccurr, "together", "solo"))

