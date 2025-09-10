packages <- c('tidyverse', 'growthrates')
lapply(packages, library, character.only=TRUE)

data <- read.csv("data/comptest/comptestdata.csv")
data$treatment <- factor(data$treatment)

