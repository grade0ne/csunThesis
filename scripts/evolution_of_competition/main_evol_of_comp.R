library(tidyverse)

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

INIT_PARAMS <- c(y0 = 1, r = 0.3, alpha = 0.2)
LOWER_PARAMS <- c(y0 = 0, r = 1e-3, alpha = 0)
UPPER_PARAMS <- c(y0 = 2, r = 8, alpha = Inf)

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
fit_models <- function(initial = INIT_PARAMS, lower = LOWER_PARAMS, upper = UPPER_PARAMS) {
  all_growthmodels(
    count ~ day | evolvedTemp + currentTemp + competition + repNum, # CUSTOMIZE TO DATA SET
    data = rPopData, # CUSTOMIZE TO DATA SET
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
  
  columns <- floor(sqrt(length(model_ids)))
  rows <- columns + 1
  par(mfrow = c(rows, columns), mar = c(1, 1, 1, 1))
  
  for (i in seq_along(models)) {
    plot(models[[i]])
    
    this_r <- round(coefs$r[i], 6)
    this_a <- round(coefs$alpha[i], 6)
    this_model <- rownames(coefs)[i]
    
    mtext(this_model, side = 3, adj = 0.1, line = -1, cex = 0.5)
    mtext(paste("r:", this_r), side = 3, adj = 0.1, line = -2, cex = 0.5)
    mtext(paste("alpha:", this_a), side = 3, adj = 0.1, line = -3, cex = 0.5)
  }
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
