lbs<-c('tidyverse','car','moments','lme4','lmerTest','growthrates','ggpattern','agricolae','emmeans','broom')
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

plotCurves <- function(model, results, rows, columns) {
  par(mfrow = c(rows, columns), mar = c(2, 2, 2, 1))
  x <- length(results$treatment)
  plotIndices <- c(1:x)
  
  for (i in plotIndices) {
    plot(model[i])
    
    this_r <- round(results$mumax[i], 4)
    this_K <- round(results$K[i], 4)
    this_temp <- results$currentTemp[i]
    this_id <- paste(results$treatment[i], results$repNum[i], sep = '-')
    
    title(main = paste(this_id, this_temp, sep = '_'), cex.main = 1)
    mtext(paste('r:', this_r), side = 3, adj = 0.15, line = -1, cex = 0.75)
    mtext(paste('K:', this_K), side = 3, adj = 0.15, line = -2, cex = 0.75)
  }
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

protistUpperP <- c(y0 = 30, mumax = 12, K = 1500)


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

rotiferResults <- rotiferResults %>% mutate(logmumax = log(mumax), dblogmumax = log(logmumax + 2))
model_rotifer_logr.ETC <- lm(logmumax ~ evolvedTemp * currentTemp * competition, rotiferResults)
model_rotifer_dblogr.ETC <- lm(dblogmumax ~ evolvedTemp * currentTemp * competition, rotiferResults)

qqp(resid(model_rotifer_logr.ETC), 'norm')
qqp(resid(model_rotifer_dblogr.ETC), 'norm')

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

emm <- emmeans(model_rotifer_logr.ETC, ~ competition | currentTemp)
pairs(emm, adjust = "tukey")

emm <- emmeans(model_rotifer_logr.ETC, ~ currentTemp| competition)
pairs(emm, adjust = "tukey")

emm <- emmeans(model_rotifer_logr.TC, ~ competition | currentTemp)
pairs(emm, adjust = "tukey")

emm <- emmeans(model_rotifer_logr.TC, ~ currentTemp| competition)
pairs(emm, adjust = "tukey")



################################################################################
# Protist r

protistResults <- protistResults %>% 
  mutate(logmumax = log(mumax),
         compFact = factor(case_when(
              evolvedTemp == '25' & competition == TRUE ~ 'rotif25',
              evolvedTemp == '30' & competition == TRUE ~ 'rotif30',
              competition == FALSE ~ 'noComp')))

model_protist_logr.CT <- lm(logmumax ~ currentTemp * compFact, protistResults)

emm <- emmeans(model_protist_logr.CT, ~ compFact | currentTemp)
pairs(emm, adjust = "tukey")

c1 <- c(-1, 0.5, 0.5)
c2 <- c(0, -1, 1)

contrastMatrix <- cbind(c1, c2)

contrasts(protistResults$compFact) <- contrastMatrix

summary(model_protist_logr.CT, split = 
          list(compFact = list("Effect of competition" = 1,
                               "Effect of rotifer evol. hist." = 2)))

# 3.87 - r no comp 30
# 2.71 - r rotif30

# 43% higher w/o compRotif30 at 30

################################################################################
# Rotifer K

model_rotifer_K.ETC <- lm(K ~ evolvedTemp * currentTemp * competition, rotiferResults)
qqp(resid(model_rotifer_K.ETC), 'norm')

model_rotifer_K.ET  <- lm(K ~ evolvedTemp * currentTemp, rotiferResults)
model_rotifer_K.EC  <- lm(K ~ evolvedTemp * competition, rotiferResults)
model_rotifer_K.TC  <- lm(K ~ currentTemp * competition, rotiferResults)
model_rotifer_K.E   <- lm(K ~ evolvedTemp, rotiferResults)
model_rotifer_K.T   <- lm(K ~ currentTemp, rotiferResults)
model_rotifer_K.C   <- lm(K ~ competition, rotiferResults)

AIC(model_rotifer_K.ETC)
AIC(model_rotifer_K.ET)
AIC(model_rotifer_K.EC) # *
AIC(model_rotifer_K.TC)
AIC(model_rotifer_K.E)
AIC(model_rotifer_K.T)
AIC(model_rotifer_K.C)


emm <- emmeans(model_rotifer_K.ETC, ~ evolvedTemp | competition | currentTemp)
pairs(emm, adjust = 'tukey')

rotiferGraph <- as.data.frame(emmeans(model_rotifer_K.ETC, ~ evolvedTemp | competition | currentTemp))

p3 <- position_dodge(width = 0.65)
ggplot(rotiferGraph, aes(x = competition, y = emmean, fill = evolvedTemp)) +
  geom_bar(stat = 'identity', position = p3, width = 0.6) +
  geom_errorbar(stat = 'identity', position = p3, width = 0.5, 
                aes(ymin = emmean - SE, ymax = emmean + SE)) +
  labs(x = 'Competition', y = 'Carrying capacity (K)', fill = 'Evol. Hist.') +
  facet_wrap(~currentTemp) +
  theme_alex

emm <- emmeans(model_rotifer_K.EC, ~ evolvedTemp | competition)
pairs(emm, adjust = 'tukey')

rotiferGraph <- as.data.frame(emmeans(model_rotifer_K.EC, ~ evolvedTemp | competition))

p3 <- position_dodge(width = 0.65)
ggplot(rotiferGraph, aes(x = competition, y = emmean, fill = evolvedTemp)) +
  geom_bar(stat = 'identity', position = p3, width = 0.6) +
  geom_errorbar(stat = 'identity', position = p3, width = 0.5, 
                aes(ymin = emmean - SE, ymax = emmean + SE)) +
  theme_alex

emm <- emmeans(model_rotifer_K.ETC, ~ evolvedTemp)
pairs(emm, adjust = 'tukey')

emm <- emmeans(model_rotifer_K.ETC, ~ competition)
pairs(emm, adjust = 'tukey')
################################################################################
# Protist K

model_protist_K.CT <- lm(K ~ currentTemp * compFact, protistResults)


emm <- emmeans(model_protist_K.CT, ~ compFact | currentTemp)
pairs(emm, adjust = "tukey")

HSD.test(model_protist_K.CT, )

protistGraph <- as.data.frame(emmeans(model_protist_K.CT, ~ compFact | currentTemp))

p3 <- position_dodge(width = 0.65)
ggplot(protistGraph, aes(x = compFact, y = emmean, fill = compFact)) +
  geom_bar(stat = 'identity', position = p3, width = 0.6) +
  geom_errorbar(stat = 'identity', position = p3, width = 0.5, 
                aes(ymin = emmean - SE, ymax = emmean + SE)) +
  geom_text(aes(label=tukey), position = p3, vjust=-2.7)+
  facet_wrap(~currentTemp) +
  theme_alex

c1 <- c(0, -1, 1)
c2 <- c(-1, 0.5, 0.5)

contrastMatrix <- cbind(c1, c2)

contrasts(protistResults$compFact) <- contrastMatrix

summary(model_protist_K.CT, split = 
          list(compFact = list("Effect of competition" = 1,
                               "Effect of rotifer evol. hist." = 2)))

contrast(emm, list(
  comp = c(-1, 0.5, 0.5),
  evol = c(0, -1, 1)
),
adjust = "tukey")

# without nocomp level

protistResults_noComp <- protistResults %>%
  filter(compFact != 'noComp')

model_protist_K.CT_NC <- lm(K ~ currentTemp * compFact, protistResults_noComp)
