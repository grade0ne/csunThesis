lbs<-c('tidyverse','car','moments','lme4','lmerTest','growthrates')
lapply(lbs,library,character.only=TRUE)
theme_alex<-readRDS('theme_alex.rds')

rawData<-read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(evolvedTemp,currentTemp,competition,species,treatID,repNum,repID),
                as.factor))

groupedData <- rawData %>% group_by(day, currentTemp, competition, treatID, species) %>% 
  summarize(mean = mean(count), ci = (sd(count)/sqrt(length(count))) * 1.96)

d1 = position_dodge(width=0.1)

################################################################################
# Protists | mean + ci x time

ggplot(groupedData%>%filter(species=='protist'), aes(x=day, y=mean, color=treatID)) +
  geom_point(position = d1) +
  geom_line(position = d1) +
  geom_errorbar(stat = 'identity', position = d1, aes(ymin = mean - ci, ymax = mean + ci)) +
  geom_point(data = rawData %>% filter(species == 'protist'), 
             aes(x = day, y = count, color = treatID), alpha = 0.15, position = d1) +
  labs(x = "Time (d)", y = "Count (indiv per 0.1 ml)") +
  theme_alex +
#  facet_wrap(~competition)
#  facet_wrap(~currentTemp)
  facet_grid(competition~currentTemp)

################################################################################
# Rotifers | mean + ci x time

ggplot(groupedData%>%filter(species=='rotifer'), aes(x=day, y=mean, color=treatID)) +
  geom_point(position = d1) +
  geom_line(position = d1) +
  geom_errorbar(stat = 'identity', position = d1, aes(ymin = mean - ci, ymax = mean + ci)) +
  geom_point(data = rawData %>% filter(species == 'rotifer'), 
             aes(x = day, y = count, color = treatID), alpha = 0.15, position = d1) +
  theme_alex +
#  facet_wrap(~competition)
#  facet_wrap(~currentTemp)
  facet_grid(competition~currentTemp)

################################################################################
# Growth Rates | Rotifers

rotiferParameters <- c(y0 = 4, mumax = 0.3, K = 100)

rotiferLowerP <- c(y0 = 0, mumax = 1e-2, K = 30)
rotiferUpperP <- c(y0 = 10, mumax = 2, K = 500)

rotiferLogistic <- all_growthmodels(
  count ~ day | treatID + evolvedTemp + currentTemp + competition + species + repNum,
  data = rawData%>%filter(species=='rotifer'),
  p = rotiferParameters,
  FUN = grow_logistic
)

rotiferResults <- as.data.frame(coef(rotiferLogistic)) %>%
  rownames_to_column(var = "id") %>%
  separate(id, into = c("treatment", "evolvedTemp", "currentTemp", "competition", "species", "repNum"), sep = ":", convert = TRUE) %>%
  mutate(treatment = factor(treatment),
         species = factor(species),
         repNum = factor(repNum),
         evolvedTemp = factor(evolvedTemp),
         currentTemp = factor(currentTemp))

rotiferGraph <- rotiferResults %>% filter(species=='rotifer') %>%
  group_by(evolvedTemp, currentTemp, competition, species) %>%
  summarize(
    mean = mean(mumax),
    ci = sd(mumax)/sqrt(length(mumax)) * 1.96
  )

d2 <- position_dodge(width=1)

ggplot(rotiferGraph, aes(x = evolvedTemp, y = mean, fill = competition)) + 
  geom_bar(stat = 'identity', position = d2) +
  geom_errorbar(stat = 'identity', position = d2, aes(ymax = mean + ci, ymin = mean - ci)) +
  facet_wrap(~currentTemp) +
  theme_alex

################################################################################
# Growth Rates | Protists

protistParameters <- c(y0 = 10, mumax = 1.5, K = 750)

protistLowerP <- c(y0 = 0, mumax = 1e-2, K = 500)
protistUpperP <- c(y0 = 50, mumax = 8, K = 1500)

protistLogistic <- all_growthmodels(
  count ~ day | treatID + currentTemp + competition + species + repNum,
  data = rawData%>%filter(species=='protist'),
  p = protistParameters,
  upper = protistUpperP,
  lower = protistLowerP,
  FUN = grow_logistic
)

protistResults <- as.data.frame(coef(protistLogistic)) %>%
  rownames_to_column(var = "id") %>%
  separate(id, into = c("treatment", "currentTemp", "competition", "species", "repNum"), sep = ":", convert = TRUE) %>%
  mutate(treatment = factor(treatment),
         species = factor(species),
         repNum = factor(repNum),
         currentTemp = factor(currentTemp))

plot_indices <- c(sample(1:60, 60))
par(mfrow = c(6, 10), mar = c(2, 2, 2, 1))

for (i in plot_indices) {
  plot(protistLogistic[[i]])
  
  this_r <- round(protistResults$mumax[i], 4)
  this_K <- round(protistResults$K[i], 4)
  this_temp <- protistResults$currentTemp[i]
  this_id <- paste(protistResults$treatment[i], protistResults$repNum[i], sep = "-")
  
  title(main = paste(this_id, this_temp, sep = "_"), cex.main = 1)
  mtext(paste("r:", this_r), side = 3, adj = 0.15, line = -1, cex = 0.75)
  mtext(paste("K:", this_K), side = 3, adj = 0.15, line = -2, cex = 0.75)
}


protistGraph <- protistResults %>%
  group_by(treatment, currentTemp, competition) %>%
  summarize(
    mean = mean(mumax),
    ci = sd(mumax)/sqrt(length(mumax)) * 1.96
  )

d2 <- position_dodge(width=1)

ggplot(protistGraph, aes(x = treatment, y = mean, fill = currentTemp)) + 
  geom_bar(stat = 'identity', position = d2) +
  geom_errorbar(stat = 'identity', position = d2, aes(ymax = mean + ci, ymin = mean - ci)) +
  facet_wrap(~competition) +
  theme_alex

ggplot(protistResults, aes(x = treatment, y = mumax, color = currentTemp)) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  theme_alex


################################################################################
# Linear Model

