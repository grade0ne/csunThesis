library(car)
library(tidyverse)


# NOV 4 ISO
data <- read.csv("data/NOV4_ISO.csv")

data$five_test <- data$count > 5

data$test_ID <- ifelse(data$five_test==TRUE, data$ID, NA)

iso_2 <- data.frame(
  ID = data$test_ID) %>%
  drop_na() %>%
  count(ID, name = "Count")

iso_2

# STOCK CHECK 4
data <- read.csv("data/stock_chk4.csv")

data$five_test <- data$count > 5

data$test_ID <- ifelse(data$five_test==TRUE, data$ID, NA)

check_4 <- data.frame(
  ID = data$test_ID) %>%
  drop_na() %>%
  count(ID, name = "Count")

check_4

# STOCK CHECK 5
data <- read.csv("data/stock_chk5.csv")

data$five_test <- data$mL > 5

data$test_ID <- ifelse(data$five_test==TRUE, data$ID, NA)

check_4 <- data.frame(
  ID = data$test_ID) %>%
  drop_na() %>%
  count(ID, name = "Count")

# STOCK CHECK 8
data <- read.csv("data/stock_chk8.csv")

data$site <- as.factor(data$site)
data$stock <- as.factor(data$stock)
data$leaf <- as.factor(data$leaf)

glimpse(data)

data$five_test <- data$countmL >= 5
data$test_ID <- ifelse(data$five_test==TRUE, data$id, NA)

check_8 <- data %>%
  filter(five_test) %>%
  group_by(id) %>%
  summarise(
    Count = n(),
    Stock_Numbers = toString(unique(stock))
  )

check_8

write.csv(check_8, "check_8_results.csv", row.names = FALSE)

