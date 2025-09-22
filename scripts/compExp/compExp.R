lbs<-c('tidyverse','car','moments')
lapply(lbs,library,character.only=TRUE)

rawData <- read.csv("data/compexp/compexpData.csv")
