library(tidyverse)

data <- read.csv("data/stocks/stock_checkup.csv") %>%
  mutate(
    stock = as.factor(stock),
    id = as.factor(id)
  ) %>%
  filter(tf==TRUE)

summary <- data %>%
  group_by(id) %>%
  summarize(count = length(id))
