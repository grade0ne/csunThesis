lbs<-c('tidyverse','car','moments','lme4','lmerTest','growthrates','ggpattern')
lapply(lbs,library,character.only=TRUE)
theme_alex<-readRDS('theme_alex.rds')

rawData<-read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(evolvedTemp,currentTemp,competition,species,treatID,repNum,repID),
                as.factor))

groupedData <- rawData %>% group_by(day, evolvedTemp, currentTemp, competition, treatID, species) %>% 
  summarize(mean = mean(count), ci = (sd(count)/sqrt(length(count))) * 1.96)

d1 = position_dodge(width=0.1)

################################################################################
# Plot of growth         | Protists  | mean + ci x time

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
# Plot of growth          | Rotifers  | mean + ci x time

ggplot(groupedData%>%filter(species=='rotifer'), aes(x=day, y=mean, color=treatID, linetype = evolvedTemp)) +
  geom_point(position = d1) +
  geom_line(position = d1) +
  scale_linetype_discrete() +
  geom_errorbar(stat = 'identity', position = d1, aes(ymin = mean - ci, ymax = mean + ci)) +
  geom_point(data = rawData %>% filter(species == 'rotifer'), 
             aes(x = day, y = count, color = treatID), alpha = 0.15, position = d1) +
  theme_alex +
#  facet_wrap(~competition)
#  facet_wrap(~currentTemp)
  facet_grid(competition~currentTemp)

################################################################################
#
#     Rotifers
#


#####################################################
# Estimating Growth Parameters | Rotifers

rotiferParameters <- c(y0 = 4, mumax = 0.3, K = 100)

rotiferLowerP <- c(y0 = 0, mumax = 1e-2, K = 30)
rotiferUpperP <- c(y0 = 10, mumax = 8, K = 500)

rotiferLogistic <- all_growthmodels(
  count ~ day | treatID + evolvedTemp + currentTemp + competition + species + repNum,
  data = rawData %>% filter(species == 'rotifer'),
  p = rotiferParameters,
  upper = rotiferUpperP,
  lower = rotiferLowerP,
  FUN = grow_logistic
)

rotiferResults <- as.data.frame(coef(rotiferLogistic)) %>%
  rownames_to_column(var = "id") %>%
  separate(id, 
           into = c("treatment", "evolvedTemp", "currentTemp", 
                        "competition", "species", "repNum"), 
           sep = ":", 
           convert = TRUE) %>%
  mutate(treatment = factor(treatment),
         species = factor(species),
         repNum = factor(repNum),
         evolvedTemp = factor(evolvedTemp),
         currentTemp = factor(currentTemp),
         competition = factor(competition))

plot_indices <- c(1:80)
par(mfrow = c(8, 10), mar = c(2, 2, 2, 1))

for (i in plot_indices) {
  plot(rotiferLogistic[[i]])
  
  this_r <- round(rotiferResults$mumax[i], 4)
  this_K <- round(rotiferResults$K[i], 4)
  this_temp <- rotiferResults$currentTemp[i]
  this_id <- paste(rotiferResults$treatment[i], rotiferResults$repNum[i], sep = "-")
  
  title(main = paste(this_id, this_temp, sep = "_"), cex.main = 1)
  mtext(paste("r:", this_r), side = 3, adj = 0.15, line = -1, cex = 0.75)
  mtext(paste("K:", this_K), side = 3, adj = 0.15, line = -2, cex = 0.75)
}


#####################################################
# Plotting Parameters | Rotifers

rotiferGraphs <- rotiferResults %>%
  group_by(evolvedTemp, currentTemp, competition, species) %>%
  summarize(
    mean_r = mean(mumax),
    ci_r = sd(mumax)/sqrt(length(mumax)) * 1.96,
    mean_K = mean(K),
    ci_K = sd(K)/sqrt(length(K)) * 1.96,
  )

d2 <- position_dodge(width=1)

# r
ggplot(rotiferGraphs, aes(x = evolvedTemp, y = mean_r, fill = competition)) + 
  geom_bar(stat = 'identity', position = d2) +
  geom_errorbar(stat = 'identity', position = d2, 
                aes(ymax = mean_r + ci_r, ymin = mean_r - ci_r)) +
  labs(x = 'Historic temperature', 
       y = 'Maximum growth rate of rotifers (r)', 
       fill = 'Competition') +
  facet_wrap(~currentTemp) +
  theme_alex

# K
ggplot(rotiferGraphs, aes(x = evolvedTemp, y = mean_K, fill = competition)) + 
  geom_bar(stat = 'identity', position = d2) +
  geom_errorbar(stat = 'identity', position = d2, 
                aes(ymax = mean_K + ci_K, ymin = mean_K - ci_K)) +
  labs(x = 'Historic temperature', 
       y = 'Carrying capacity of rotifers (K)', 
       fill = 'Competition') +
  facet_wrap(~currentTemp) +
  theme_alex

################################################################################
#
#     Protists
#


#####################################################
# Estimating Growth Parameters | Protists


protistParameters <- c(y0 = 10, mumax = 1.5, K = 750)

protistLowerP <- c(y0 = 0, mumax = 1e-2, K = 500)
protistUpperP <- c(y0 = 30, mumax = 8, K = 1500)

protistLogistic <- all_growthmodels(
  count ~ day | treatID + evolvedTemp + currentTemp + competition + species + repNum,
  data = rawData %>% filter(species == 'protist'),
  p = protistParameters,
  upper = protistUpperP,
  lower = protistLowerP,
  FUN = grow_logistic
)

protistResults <- as.data.frame(coef(protistLogistic)) %>%
  rownames_to_column(var = "id") %>%
  separate(id, 
           into = c("treatment", "evolvedTemp", "currentTemp", "competition", 
                        "species", "repNum"), 
           sep = ":", 
           convert = TRUE) %>%
  mutate(treatment = factor(treatment),
         species = factor(species),
         repNum = factor(repNum),
         currentTemp = factor(currentTemp),
         evoRotifTemp = case_when(
           treatment %in% c('HP', 'NP') ~ "Protist only (25C)",
           TRUE ~ paste0(evolvedTemp, "C")
         ))

plot_indices <- c(1:60)
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

#####################################################
# Plotting Parameters | Protists

protistGraphs <- protistResults %>%
  group_by(treatment, evoRotifTemp, currentTemp, competition) %>%
  summarize(
    mean_r = mean(mumax),
    ci_r = sd(mumax) / sqrt(length(mumax)) * 1.96,
    mean_K = mean(K),
    ci_K = sd(K) / sqrt(length(K)) * 1.96
  )

d2 <- position_dodge(width=1)

# r
ggplot(protistGraphs, aes(x = treatment, y = mean_r, fill = currentTemp, color = evoRotifTemp)) + 
  geom_bar(stat = 'identity', position = d2, linewidth = 1) +
  geom_errorbar(stat = 'identity', position = d2, aes(ymax = mean_r + ci_r, ymin = mean_r - ci_r)) +
  labs(x = 'Treatment', y = 'Maximum growth rate of protists (r)', fill = 'Contemp. Temp') +
  scale_fill_manual(values = c('lightgreen', 'orange')) +
  scale_color_manual(values = c("cyan3", "red2", "black")) +
  theme_alex

# k
ggplot(protistGraphs, aes(x = treatment, y = mean_K, fill = currentTemp, color = evoRotifTemp)) + 
  geom_bar(stat = 'identity', position = d2, linewidth = 1.4) +
  geom_errorbar(stat = 'identity', position = d2, aes(ymax = mean_K + ci_K, ymin = mean_K - ci_K)) +
  labs(x = 'Treatment', y = 'Carrying capacity of protists (K)', fill = 'Temperature') +
  scale_fill_manual(values = c('#90C3F2', '#FFB000')) +
  scale_color_manual(values = c("#648FFF", "#DE7604", "black")) +
  theme_alex

################################################################################
# ANOVA
#
#   Dvs: rotifer r, protist r, rotifer K, protist K, size <- each own 3-way ANOVA
#
#   Evol hist
#   Cont temp
#   Comp
#

################################################################################
# Rotifer r

model_rotifer_r.ETC <- lm(mumax ~ evolvedTemp * currentTemp * competition, rotiferResults)
qqp(resid(model_rotifer_r.ETC), 'norm')

rotiferResults <- rotiferResults %>% mutate(logmumax = log(mumax))
model_rotifer_logr.ETC <- lm(logmumax ~ evolvedTemp * currentTemp * competition, rotiferResults)
qqp(resid(model_rotifer_logr.ETC), 'norm')

model_rotifer_logr.ET  <- lm(logmumax ~ evolvedTemp * currentTemp, rotiferResults)
model_rotifer_logr.EC  <- lm(logmumax ~ evolvedTemp * competition, rotiferResults)
model_rotifer_logr.TC  <- lm(logmumax ~ currentTemp * competition, rotiferResults)
model_rotifer_logr.E   <- lm(logmumax ~ evolvedTemp, rotiferResults)
model_rotifer_logr.T   <- lm(logmumax ~ currentTemp, rotiferResults)
model_rotifer_logr.C   <- lm(logmumax ~ competition, rotiferResults)


AIC(model_rotifer_logr.ETC) # 
AIC(model_rotifer_logr.ET)  # 
AIC(model_rotifer_logr.EC)  # 
AIC(model_rotifer_logr.TC)  #  *
AIC(model_rotifer_logr.E)   # 
AIC(model_rotifer_logr.T)   # 
AIC(model_rotifer_logr.C)   # 

################################################################################
# Protist r

model_protist_r.ETC <- lm(mumax ~ evolvedTemp * currentTemp * competition, protistResults)
qqp(resid(model_protist_r.ETC), 'norm')

protistResults <- protistResults %>% mutate(logmumax = log(mumax))

model_protist_logr.ETC <- lm(logmumax ~ evolvedTemp * currentTemp * competition, protistResults)
qqp(resid(model_protist_logr.ETC), 'norm')

model_protist_logr.ET  <- lm(logmumax ~ evolvedTemp * currentTemp, protistResults)
model_protist_logr.EC  <- lm(logmumax ~ evolvedTemp * competition, protistResults)
model_protist_logr.TC  <- lm(logmumax ~ currentTemp * competition, protistResults)
model_protist_logr.E   <- lm(logmumax ~ evolvedTemp, protistResults)
model_protist_logr.T   <- lm(logmumax ~ currentTemp, protistResults)
model_protist_logr.C   <- lm(logmumax ~ competition, protistResults)


AIC(model_protist_logr.ETC) # 
AIC(model_protist_logr.ET)  # 
AIC(model_protist_logr.EC)  # 
AIC(model_protist_logr.TC)  #  *
AIC(model_protist_logr.E)   # 
AIC(model_protist_logr.T)   # 
AIC(model_protist_logr.C)   # 

################################################################################
# Rotifer K



################################################################################
# Protist K

model_protist_K.ETC <- lm(K ~ evolvedTemp * currentTemp * competition, protistResults)
qqp(resid(model_protist_K.ETC), 'norm')

model_protist_K.ET  <- lm(K ~ evolvedTemp * currentTemp, protistResults)
model_protist_K.EC  <- lm(K ~ evolvedTemp * competition, protistResults)
model_protist_K.TC  <- lm(K ~ currentTemp * competition, protistResults)
model_protist_K.E   <- lm(K ~ evolvedTemp, protistResults)
model_protist_K.T   <- lm(K ~ currentTemp, protistResults)
model_protist_K.C   <- lm(K ~ competition, protistResults)


AIC(model_protist_K.ETC) # 
AIC(model_protist_K.ET)  # 
AIC(model_protist_K.EC)  # 
AIC(model_protist_K.TC)  #   *
AIC(model_protist_K.E)   # 
AIC(model_protist_K.T)   # 
AIC(model_protist_K.C)   # 
