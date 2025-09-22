lbs<-c('tidyverse','car','moments','lme4','lmerTest')
lapply(lbs,library,character.only=TRUE)
rawData<-read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(replicate,evolvedTemp,currentTemp,competition,species,idCode),
                as.factor))