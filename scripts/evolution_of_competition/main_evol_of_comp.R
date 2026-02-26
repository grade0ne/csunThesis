library(tidyverse)

#### DATA HANDLING ####

popData <- read.csv("data/compexp/compexpData.csv") %>%
  mutate(hour = as.integer(day * 24),
         across(c(evolvedTemp, currentTemp, competition, repNum), factor))

rPopData <- popData %>%
  filter(species == "rotifer")

pPopData <- popData %>%
  filter(species == "protist")


#### GROWTH MODELS ####

library(growthrates)

# Rotifer parameters: (TIME IN HOURS)
INIT_PARAMS <- c(y0 = 5, r = 0.1, alpha = 0.0005)
LOWER_PARAMS <- c(y0 = 3, r = 0.003, alpha = 0.0001)
UPPER_PARAMS <- c(y0 = 7, r = 8, alpha = 1)

# Protist parameters:
INIT_PARAMS <- c(y0 = 10, r = 0.3, alpha = 0.15)
LOWER_PARAMS <- c(y0 = 0, r = 0.003, alpha = 0.001)
UPPER_PARAMS <- c(y0 = 80, r = 20, alpha = 1)

# Protist parameters: (TIME IN HOURS)
INIT_PARAMS <- c(y0 = 10, r = 0.1, alpha = 0.0003)
LOWER_PARAMS <- c(y0 = 5, r = 0.001, alpha = 0.000005)
UPPER_PARAMS <- c(y0 = 15, r = 5, alpha = 0.02)

# Defining custom function per Part 2, sec. 4 of growthrates documentation
logistic_alpha_function <- function(time, parms) {
  y0    <- parms[["y0"]]
  r     <- parms[["r"]]
  alpha <- parms[["alpha"]]
  
  y <- (r * y0) / (alpha * y0 + (r - alpha * y0) * exp(-r * time))
  
  cbind(time = time, y = y)
}

grow_logistic_alpha <- growthmodel(
  logistic_alpha_function,
  pnames = c("y0", "r", "alpha")
)


# Fit logistic growth model to growth data set. Returns multiple nonlinear fits.
# usage: 'model_object_name' <- fit_models()
fit_models <- function(data, initial = INIT_PARAMS, lower = LOWER_PARAMS, upper = UPPER_PARAMS) {
  all_growthmodels(
    count ~ hour | evolvedTemp + currentTemp + competition + repNum, # CUSTOMIZE TO DATA SET
    data = data, # CUSTOMIZE TO DATA SET
    p = initial,
    upper = upper,
    lower = lower,
    FUN = grow_logistic_alpha
  )
}

# Visually inspect fit of models
plot_models <- function(models) {
  coefs <- as.data.frame(coef(models))
  model_ids <- c(1:length(coefs$y0))
  
  columns <- ceiling(sqrt(length(model_ids)))
  rows <- columns + 1
  par(mfrow = c(rows, columns), mar = c(1, 1, 1, 1))
  
  for (i in seq_along(model_ids)) {
    
    plot(models[[i]])
    
    this_r <- round(coefs$r[i], 6)
    this_a <- round(coefs$alpha[i], 6)
    this_model <- rownames(coefs)[i]
    
    mtext(this_model, side = 3, adj = 0.1, line = -1, cex = 0.5)
    mtext(paste("r:", this_r), side = 3, adj = 0.1, line = -2, cex = 0.5)
    mtext(paste("alpha:", this_a), side = 3, adj = 0.1, line = -3, cex = 0.5)
  }
  
  par(mfrow = c(1, 1))
}


#### ANOVA ####

rGrowData <- results(r_models) %>%
  rownames_to_column(var = "id") %>%
  separate(col = "id", into = c("evolvedTemp", "currentTemp", "competition", "replicate"), sep = ":") %>%
  mutate(across(c(evolvedTemp, currentTemp, competition, replicate), factor)) %>%
  mutate(r_hour = r,
         r_day = r * 24,
         alpha_hour = alpha,
         alpha_day = alpha * 24,
         K = r_day / alpha_day)

pGrowData <- results(p_models) %>%
  rownames_to_column(var = "id") %>%
  separate(col = "id", into = c("evolvedTemp", "currentTemp", "competition", "replicate"), sep = ":") %>%
  mutate(across(c(evolvedTemp, currentTemp, competition, replicate), factor),
         compFactor = case_when(evolvedTemp == "25" & competition == TRUE ~ "25-evol-rotif",
                                evolvedTemp == "30" & competition == TRUE ~ "30-evol-rotif",
                                competition == FALSE ~ "no-comp")) %>%
  mutate(r_hour = r,
         r_day = r * 24,
         alpha_hour = alpha,
         alpha_day = alpha * 24)

library(car)
library(MASS)
library(moments)
library(emmeans)

# Rotifer r
model1 <- lm(log(r_day) ~ evolvedTemp * currentTemp * competition, data = rGrowData)
Anova(model1, type = "III")

emm1 <- emmeans(model1, ~ evolvedTemp | competition)
pairs(emm1, adjust = "tukey")


# Rotifer a
model2 <- lm(log(alpha_day) ~ evolvedTemp * currentTemp * competition, data = rGrowData)
Anova(model2, type = "III")

emm2 <- emmeans(model2, ~ competition | evolvedTemp)
pairs(emm2, adjust = "tukey")


# Protist r
model3 <- lm(log(r_day) ~ currentTemp * compFactor, data = pGrowData)
Anova(model3, type = "III")

emm3 <- emmeans(model3, ~ compFactor | currentTemp)
pairs(emm3, adjust = "tukey")


# Protist a
model4 <- lm(log(alpha_day) ~ currentTemp * compFactor, data = pGrowData)
Anova(model4, type = "III")

emm4 <- emmeans(model4, ~ compFactor | currentTemp)
pairs(emm4, adjust = "tukey")


#### FIGURES ####

fig1data <- rGrowData %>%
  group_by(evolvedTemp, currentTemp, competition) %>%
  summarize(mean_r = mean(r_day), ci_r = 1.96 * sd(r_day)/sqrt(length(r_day)))

ggplot(fig1data, aes(x = currentTemp, y = mean_r, fill = competition)) +
  geom_bar(stat = "identity", position= position_dodge(width = 1)) +
  geom_errorbar(stat = "identity", position = position_dodge(width = 1),
                aes(ymin = mean_r - ci_r, ymax = mean_r + ci_r)) +
  facet_wrap(~evolvedTemp) +
  theme_minimal()





