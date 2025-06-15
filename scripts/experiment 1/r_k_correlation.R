library(car)
library(tidyverse)
library(dplyr)
library(moments)

data <- read.csv("data/log_results.csv")

plot(K~mumax, data)

data2 <- data %>%
  filter(K < 199, K > 30.01) %>%
  filter(mumax < 0.599) %>%
  mutate(
    logR = log(mumax),
    logK = log(K)
  )

plot(K~mumax, data2)

qqp(data2$mumax, "norm")
qqp(data2$K, "norm")

kurtosis(data2$mumax)
kurtosis(data2$K)

plot(logK~logR, data2)

qqp(data2$logR, "norm")
qqp(data2$logK, "norm")

cor.test(data2$logR, data2$logK, method="pearson", na.rm=TRUE)
cor.test(data2$logR, data2$logK, method="spearman", na.rm=TRUE)


#########################
# By clone

data2$clone <- as.factor(data2$clone)

data3 <- data2 %>%
  group_by(clone) %>%
  summarise(
    mean_logK = mean(logK),
    mean_logR = mean(logR)
  )

plot(data3$mean_logK~data3$mean_logR)

cor.test(data3$mean_logK, data3$mean_logR, method = "pearson")

#########################
# Multiple regression

model1 <- lm(logK ~ logR * treatment, data = data2)

plot(model1)
anova(model1)


predicted <- expand.grid(
  logR = seq(min(data2$logR), max(data2$logR), length.out = 100),
  treatment = unique(data2$treatment)
)

predicted$predicted_logK <- predict(model1, newdata = predicted)

ggplot(data2, aes(x = logR, y = logK, color = treatment)) +
  geom_point() +
  geom_line(data = predicted, aes(x = logR, y = predicted_logK, color = treatment)) +
  theme_minimal()

#########################
# Plots

ggplot(data2, aes(x = logR, y = logK, color = treatment)) +
  geom_point(stat = "identity")

grouped_data <- data2 %>%
  group_by(treatment, clone) %>%
  mutate(
    meanK = mean(logK),
    meanR = mean(logR),
    seK = sd(logK)/sqrt(n()),
    seR = sd(logR)/sqrt(n())
  )


