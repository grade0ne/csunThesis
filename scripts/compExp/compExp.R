lbs<-c('tidyverse','car','moments','lme4','lmerTest')
lapply(lbs,library,character.only=TRUE)
theme_alex<-readRDS('theme_alex.rds')

rawData<-read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(evolvedTemp,currentTemp,competition,species,treatID,repNum,repID),
                as.factor))

groupedData <- rawData %>% group_by(day, currentTemp, competition, treatID, species) %>% 
  summarize(mean = mean(count), se = sd(count)/sqrt(length(count)))

d1 = position_dodge(width=0.1)

################################################################################
# Protists | mean + SE x time

ggplot(groupedData%>%filter(species=='protist'), aes(x=day, y=mean, color=treatID)) +
  geom_point(position = d1) +
  geom_line(position = d1) +
  geom_errorbar(stat = 'identity', position = d1, aes(ymin = mean - se, ymax = mean + se)) +
  geom_point(data = rawData %>% filter(species == 'protist'), 
             aes(x = day, y = count, color = treatID), alpha = 0.15, position = d1) +
  theme_alex #+
#  facet_wrap(~competition)
#  facet_wrap(~currentTemp)

################################################################################
# Rotifers | mean + SE x time

ggplot(groupedData%>%filter(species=='rotifer'), aes(x=day, y=mean, color=treatID)) +
  geom_point(position = d1) +
  geom_line(position = d1) +
  geom_errorbar(stat = 'identity', position = d1, aes(ymin = mean - se, ymax = mean + se)) +
  geom_point(data = rawData %>% filter(species == 'rotifer'), 
             aes(x = day, y = count, color = treatID), alpha = 0.15, position = d1) +
  theme_alex #+
#  facet_wrap(~competition)
#  facet_wrap(~currentTemp)
