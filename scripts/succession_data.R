library(car)
library(tidyverse)

data <- read.csv("data/PP_succession_censusdata.csv")


rotifers_vol <- data %>%
  filter(Rotifers > 0) %>%
  mutate(
    Volume..mL. = as.numeric(Volume..mL.),
    Rotifers = as.numeric(Rotifers)
  ) %>%
  filter(!is.na(Volume..mL.), !is.na(Rotifers)) %>%
  summarise(
    meanVol = mean(Volume..mL.),
    se = sd(Volume..mL.)/sqrt(n())
  ) %>%
  mutate(group = "Rotifers")

all_vol <- data %>%
  mutate(
    Volume..mL. = as.numeric(Volume..mL.)
  ) %>%
  filter(!is.na(Volume..mL.)) %>%
  summarise(
    meanVol = mean(Volume..mL.),
    se = sd(Volume..mL.)/sqrt(n())
  ) %>%
  mutate(group = "All")

merged_vol <- bind_rows(rotifers_vol, all_vol)

ggplot(merged_vol, aes(x = group, y = meanVol)) +
  geom_bar(stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin=meanVol-se, ymax=meanVol+se))+
  theme_minimal()

leaf_age <- data %>%
  filter(Rotifers > 0) %>%
  mutate(
    Plant.ID.Numeric = as.factor(Plant.ID.Numeric)
  ) %>%
  group_by(Plant.ID.Numeric) %>%
  summarise(maxAge = max(leaf_age.days.))

# switch from doing plant ID numeric to band number!

max_vol <- data %>%
  filter(!is.na(Volume..mL.)) %>%
  mutate(Band.. = as.factor(Band..),
         Volume..mL. = as.numeric(Volume..mL.)) %>%
  group_by(Band..) %>%
  summarise(
    maxVol = max(Volume..mL., na.rm = TRUE)
  )
