library(MoMAColors)
library(MetBrewer)
theme_alex<-readRDS('theme_alex.rds')

#ffb651 yellow
#2e92a2 blue
#aa3f5d berg
#61bea4 teal

# df objects and library depdenancies from scripts/compExp.R

#### Table 1 #### 
rotiferR <- tidy(Anova(model_rotifer_logr.ETC, type = 'III')) %>%
  filter(term != '(Intercept)') %>%
  mutate(across(-term, ~ round(., 4)))
write.table(rotiferR, file = "rotiferR.txt", sep = ',', quote = F, row.names = F)


#### Table 2a #### 
protistR <- tidy(Anova(model_protist_logr.CT, type = 'III')) %>%
  filter(term != '(Intercept)') %>%
  mutate(across(-term, ~ round(., 4)))
write.table(protistR, file = "protistR.txt", sep = ',', quote = F, row.names = F)


#### Table 2b #### 
emm <- emmeans(model_protist_logr.CT, ~ compFact | currentTemp)
protistRComparisons <- tidy(pairs(emm, adjust = "tukey")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.table(protistRComparisons, file = "protistRcomparisons.txt", sep = ',', quote = F, row.names = F)

#### Figure 1 #### 
fig1Data <- as.data.frame(emmeans(model_rotifer_logr.ETC, 
                                  ~ evolvedTemp)) %>%
  mutate(
    btMean = exp(emmean),
    btLowerCL = exp(lower.CL),
    btUpperCL = exp(upper.CL))


d1 <- position_dodge(width = .75)

ggplot(fig1Data, aes(x = evolvedTemp, y = btMean, fill = evolvedTemp)) +
  geom_bar(stat = 'identity', 
           position = d1, 
           width = 0.6) +
  geom_errorbar(stat = 'identity', 
                position = d1, 
                aes(ymin = btLowerCL, ymax = btUpperCL),
                width = 0.5) +
  labs(x = 'Evolved Temperature (°C)', y = 'Maximum growth rate (r ± SEM)') +
  scale_fill_manual(values = c('#2e92a2', '#ffb651')) +
  theme_alex +
  theme(legend.position = 'none')

#### Figure 2 #### 
fig2Data <- as.data.frame(emmeans(model_rotifer_logr.ETC, 
                                  ~ currentTemp * competition)) %>%
  mutate(
    btMean = exp(emmean),
    btLowerCL = exp(lower.CL),
    btUpperCL = exp(upper.CL))

d2 <- position_dodge(width = .2)

ggplot(fig2Data, aes(x = currentTemp, y = btMean, color = competition, group = competition)) +
  geom_point(stat = 'identity', 
             position = d2,
             size = 5) +
  geom_line(stat = 'identity', position = d2, linewidth = 1.2) +
  geom_errorbar(stat = 'identity', 
                position = d2, 
                aes(ymin = btLowerCL, ymax = btUpperCL),
                width = 0.3) +
  labs(x = 'Current Temperature (°C)', y = 'Maximum growth rate (r ± SEM)') +
  scale_color_manual(values = c('#61bea4', '#aa3f5d')) +
  theme_alex +
  theme(legend.position = 'none')

#### Figure 3 ####
fig3Data <- as.data.frame(emmeans(model_protist_logr.CT, 
                                  ~ currentTemp * compFact)) %>%
  mutate(
    btMean = exp(emmean),
    btLowerCL = exp(lower.CL),
    btUpperCL = exp(upper.CL))

d3 <- position_dodge(width = .15)

ggplot(fig3Data, aes(x = currentTemp, y = btMean, color = compFact, group = compFact)) +
  geom_point(stat = 'identity', 
             position = d3,
             size = 5) +
  geom_line(stat = 'identity', position = d3, linewidth = 1.2) +
  geom_errorbar(stat = 'identity', 
                position = d3, 
                aes(ymin = btLowerCL, ymax = btUpperCL),
                width = 0.3) +
  labs(x = 'Current Temperature (°C)', y = 'Maximum growth rate (r ± SEM)') +
  scale_color_manual(values = c('grey', '#2e92a2', '#ffb651')) +
  theme_alex

ggplot(fig3Data, aes(x = currentTemp, y = btMean, fill = compFact)) +
  geom_bar(stat = 'identity', 
             position = d1,
             width = 0.5) +
  geom_errorbar(stat = 'identity', 
                position = d1, 
                aes(ymin = btLowerCL, ymax = btUpperCL),
                width = 0.4) +
  labs(x = 'Current Temperature (°C)', y = 'Maximum growth rate (r ± SEM)') +
  scale_fill_manual(values = c('grey', '#2e92a2', '#ffb651')) +
  theme_alex
