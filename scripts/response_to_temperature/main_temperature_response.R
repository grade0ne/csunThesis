library(tidyverse)

#### DATA HANDLING ####

popData <- read.csv("Data/tempres/tempRes2-30d.csv") %>%
  mutate(Day = Day - 3,
         across(c(Site, Leaf, Clone, Treatment, ID), as.factor),
         across(c(Day, Count), as.integer))

clone_key <- c(
  "1" = 1,  "2" = 2,  "17" = 3,  "18" = 4,  "19" = 5,
  "72" = 6, "73" = 7, "74" = 8, "39" = 9,  "40" = 10,
  "42" = 11,"75" = 12,"26" = 13,"25" = 14,"35" = 15,
  "36" = 16,"37" = 17
)

sizeData <- read.csv("Data/sizeData/diversitySizeData.csv") %>%
  separate(base, into = c("unk", "culture", "image", "object"), sep = "_", convert = TRUE) %>%
  separate(culture, into = c("siteLeaf", "clone", "replicate"), sep = "-", convert = TRUE) %>%
  mutate(site = as.integer(case_when(str_sub(siteLeaf, 1, 1) == "c" ~ 1,
                                     str_sub(siteLeaf, 1, 1) == "p" ~ 2)),
         leaf = str_sub(siteLeaf, start = 2),
         leaf = as.integer(case_when(leaf == "23" ~ 1,
                                     leaf == "7"  ~ 2,
                                     leaf == "36" & site == "1" ~ 3,
                                     leaf == "36" & site == "2" ~ 4,
                                     leaf == "24" ~ 5,
                                     leaf == "3"  ~ 6)),
         clone = as.integer(unname(clone_key[as.character(clone)])),
         temperature = as.factor(case_when(replicate %in% c(1, 2, 3) ~ "25C",
                                           replicate %in% c(4, 5, 6) ~ "30C"))
         ) %>%
  select(
    temperature, site, leaf, clone, area_px
)

#### GROWTH MODELS ####

library(growthrates)

INIT_PARAMS <- c(y0 = 1, r = 0.3, alpha = 0.2)
LOWER_PARAMS <- c(y0 = 0, r = 0.001, alpha = 0.001)
UPPER_PARAMS <- c(y0 = 5, r = 8, alpha = 1)

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

# Test data for custom alpha function
# it's fun to play with alpha--use for defense explanation of logistic?
time <- 1:30
out <- grow_logistic_alpha(time, parms = list(y0 = 1, r = 0.5, alpha = 0.5))
plot(out[, "time"], out[, "y"], type = "b")



# Fit logistic growth model to growth data set. Returns multiple nonlinear fits.
# usage: 'model_object_name' <- fit_models()
fit_models <- function(initial = INIT_PARAMS, lower = LOWER_PARAMS, upper = UPPER_PARAMS) {
  all_growthmodels(
    Count ~ Day | Treatment + Site + Leaf + Clone + ID, # CUSTOMIZE TO DATA SET
    data = popData, # CUSTOMIZE TO DATA SET
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
  par(mfrow = c(columns, rows), mar = c(1, 1, 1, 1))
  
  for (i in seq_along(model_ids)) {
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
