library(tidyverse)

data <- read.csv("data/evol_select.csv")

data <- data %>%
  group_by(Source) %>%
  mutate(
     Cs = Count * 10,
     Vs_.25K = 437 / Cs,
     V10 = Vs * 10,
     
  )
