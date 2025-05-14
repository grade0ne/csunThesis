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

# using subset data
subset <- read.csv("data/sp_subset.csv")

subset <- subset %>%
  mutate(Rotifers = as.logical(Rotifers))

max_age <- subset %>%
  group_by(Leaf) %>%
  summarise(max_age = max(Age))

max_vol <- subset %>%
  filter(Volume > 0) %>%
  group_by(Leaf) %>%
  slice_max(order_by = Volume, n = 1, with_ties = F) %>%
  left_join(max_age, by = "Leaf") %>%
  mutate(age_perc = (Age / max_age) * 100) %>%
  select(Leaf, Age, Volume, max_age, age_perc) %>%
  mutate(group = "All")


max_vol_r <- subset %>%
  filter(Volume > 0, Rotifers == TRUE) %>%
  group_by(Leaf) %>%
  slice_max(order_by = Volume, n = 1, with_ties = F) %>%
  left_join(max_age, by = "Leaf") %>%
  mutate(age_perc = (Age / max_age) * 100) %>%
  select(Leaf, Age, Volume, max_age, age_perc) %>%
  mutate(group = "Rotifers")

vol_graph <- bind_rows(max_vol, max_vol_r)
vol_graph <- vol_graph %>%
  group_by(group) %>%
  summarise(
    Volume = Volume,
    meanVol = mean(Volume),
    seV = sd(Volume)/sqrt(n()),
    Age = mean(Age),
    seA = sd(Age)/sqrt(n())
  )

ggplot(vol_graph, aes(x = group, y = meanVol)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = meanVol - seV,  ymax = meanVol + seV), position = position_dodge(0.5)) +
  theme_minimal()
