rData <- read.csv("Data/compexp/rotiferParameters.csv")
pData <- read.csv("Data/compexp/protistParameters.csv")

rData <- rData %>%
  mutate(alpha = mumax / K,
         r = mumax,
         lnr = log(r),
         across(c(treatment, evolvedTemp, currentTemp, competition), as.factor))

plot(r~K, rData)
plot(r~alpha, rData)

plot(alpha~K, rData)

#evo

evoAlpha <- rData %>%
  group_by(evolvedTemp) %>%
  summarise(mean = mean(alpha),
            se = sd(alpha)/sqrt(length(alpha)))

evor <- rData %>%
  group_by(evolvedTemp) %>%
  summarise(mean = mean(r),
            se = sd(r)/sqrt(length(r)))

evoK <- rData %>%
  group_by(evolvedTemp) %>%
  summarise(mean = mean(K),
            se = sd(K)/sqrt(length(K)))

ggplot(evoAlpha, aes(evolvedTemp, mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(stat='identity', aes(ymin=mean-se, ymax=mean+se))

ggplot(evor, aes(evolvedTemp, mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(stat='identity', aes(ymin=mean-se, ymax=mean+se))

ggplot(evoK, aes(evolvedTemp, mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(stat='identity', aes(ymin=mean-se, ymax=mean+se))



#temp

tempAlpha <- rData%>%
  group_by(currentTemp) %>%
  summarise(mean = mean(alpha),
            se = sd(alpha)/sqrt(length(alpha)))

tempr <- rData %>%
  group_by(currentTemp) %>%
  summarise(mean = mean(r),
            se = sd(r)/sqrt(length(r)))

tempK <- rData %>%
  group_by(currentTemp) %>%
  summarise(mean = mean(K),
            se = sd(K)/sqrt(length(K)))

ggplot(tempAlpha, aes(currentTemp, mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(stat='identity', aes(ymin=mean-se, ymax=mean+se))

ggplot(tempr, aes(currentTemp, mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(stat='identity', aes(ymin=mean-se, ymax=mean+se))

ggplot(tempK, aes(currentTemp, mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(stat='identity', aes(ymin=mean-se, ymax=mean+se))


# 
model1 <- lm(log(alpha) ~ evolvedTemp * currentTemp * competition, rData)

plot(model1)
library(car)
qqp(model1)

plot(resid(model1)~fitted(model1))
abline(h=0)

Anova(model1, type="III")

#
model2 <- lm(log(r) ~ evolvedTemp * currentTemp * competition, rData)

plot(model2)
Anova(model2, type="III")
qqp(model2)

Anova(model2, type="III")

compTempR <- rData %>%
  group_by(currentTemp, competition) %>%
  summarise(mean = mean(lnr),
            se = sd(lnr)/sqrt(length(lnr)))

ggplot(compTempR, aes(currentTemp, mean, fill = competition)) +
  geom_bar(stat='identity', position='dodge') +
  geom_errorbar(stat='identity', position='dodge', aes(ymin=mean-se, ymax=mean+se)) +
  theme_minimal()

#
model3 <- lm(K ~ evolvedTemp * currentTemp * competition, rData)

plot(model3)
qqp(model3, 'norm')
Anova(model3, type="III")
