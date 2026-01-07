# Created by Alex Mendelson on Jan 7 2026

library(tidyverse)
library(lme4)
library(car)
library(lmerTest)
library(moments)

rawData <- read.csv("Data/sizeData/diversitySizeData.csv")

data <- rawData %>%
  separate(base,
           into = c("index", "id", "frame", "individual"),
           sep = "_",
           convert = TRUE) %>%
  separate(id,
           into = c("siteLeaf", "clone", "treatment_class"),
           sep = "-",
           convert = TRUE) %>%
  mutate(temperature = case_when(treatment_class %in% c("1", "2", "3") ~ 25, TRUE ~ 30)) %>%
  mutate(site = str_extract(as.character(siteLeaf), "^[A-Za-z]"),
         leaf = dense_rank(as.character(siteLeaf))) %>%
  group_by(clone) %>%
  mutate(clone_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(across(c(site, leaf, clone_id, temperature), as.factor)) %>%
  select(site, leaf, clone_id, temperature, area_px)

summary <- data %>%
  group_by(clone_id, temperature) %>%
  summarise(mean = mean(area_px),
            n = length(area_px),
            se = sd(area_px)/sqrt(length(area_px)))

plot(area_px~site, data)
plot(area_px~leaf, data)
plot(area_px~clone_id, data)
plot(area_px~temperature, data)

model1 <- lmer(area_px~temperature + (1|site/leaf/clone_id), data)

model2 <- lmer(log(area_px)~temperature + (1|site/leaf/clone_id), data)
model3 <- lm(log(area_px)~temperature, data)
anova(model2, model3)
