library(tidyverse)

#### DATA HANDLING ####

popData <- read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(day, count), as.integer),
         across(c(evolvedTemp, currentTemp, competition, repNum), factor)) %>%
  select(day, species, evolvedTemp, currentTemp, competition, repNum, count)

rPopData <- popData %>%
  filter(species == "rotifer")

pPopData <- popData %>%
  filter(species == "protist")


#### GROWTH MODELS ####

library(growthrates)

# Rotifer parameters:
INIT_PARAMS <- c(y0 = 1, r = 0.3, alpha = 0.15)
LOWER_PARAMS <- c(y0 = 0, r = 0.003, alpha = 0.003)
UPPER_PARAMS <- c(y0 = 10, r = 8, alpha = 1)

# Protist parameters:
INIT_PARAMS <- c(y0 = 10, r = 0.3, alpha = 0.15)
LOWER_PARAMS <- c(y0 = 0, r = 0.003, alpha = 0.001)
UPPER_PARAMS <- c(y0 = 80, r = 10, alpha = 1)

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
    count ~ day | evolvedTemp + currentTemp + competition + repNum, # CUSTOMIZE TO DATA SET
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

# model summaries
summarize_models <- function(models) {
  
  ids <- rownames(coef(models))
  out_list <- list()
  
  for (i in seq_along(ids)) {
    fit_id <- ids[i]
    fit <- models[[fit_id]]
    
    summary <- summary(fit)
    table <- as.data.frame(summary$par)
    
    table$term <- rownames(table)
    table$model_id <- fit_id
    rownames(table) <- NULL
    
    out_list[[i]] <- table
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

#### ANOVA ####

rGrowData <- results(models) %>%
  rownames_to_column(var = "id") %>%
  separate(col = "id", into = c("evolvedTemp", "currentTemp", "competition", "replicate"), sep = ":") %>%
  mutate(across(c(evolvedTemp, currentTemp, competition, replicate), factor)) %>%
  select(!repNum)

pGrowData <- results(models) %>%
  rownames_to_column(var = "id") %>%
  separate(col = "id", into = c("evolvedTemp", "currentTemp", "competition", "replicate"), sep = ":") %>%
  mutate(across(c(evolvedTemp, currentTemp, competition, replicate), factor),
         compFactor = case_when(evolvedTemp == "25" & competition == TRUE ~ "25-evol-rotif",
                                evolvedTemp == "30" & competition == TRUE ~ "30-evol-rotif",
                                competition == FALSE ~ "no-comp"))

library(car)
library(MASS)
library(moments)

model1 <- lm(log(r) ~ evolvedTemp * currentTemp * competition, data = rGrowData)

model2 <- lm(log(alpha) ~ evolvedTemp * currentTemp * competition, data = rGrowData)

model3 <- lm(log(r) ~ currentTemp * compFactor, data = pGrowData)

model4 <- lm(log(alpha) ~ currentTemp * compFactor, data = pGrowData)
#### FIGURES ####


