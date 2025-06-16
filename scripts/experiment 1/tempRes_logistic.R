# Growthrates package documentation
# https://cran.r-project.org/web/packages/growthrates/vignettes/Introduction.html

# Worked example of triple nested model structure - Mick Keough
# https://mjkeough.github.io/examples/caballes.nb.html#

library(afex)
library(car)
library(lattice)
library(lme4)
library(lmerTest)
library(nlme)
library(VCA)
library(Rmisc)
library(growthrates)
library(tidyverse)
library(moments)
library(MASS)

#######################################################
# Data

my_data <- read.csv("data/tempRes2-30d.csv")

my_data <- my_data %>%
  mutate(
  Count1 = my_data$Count + 1,
  Replicate = as.integer(substr(my_data$ID, nchar(my_data$ID), nchar(my_data$ID))),
  log_count = log(Count1)
)

my_data[c("Site", "Leaf", "Clone", "Treatment", "Replicate")] <-
  lapply(my_data[c("Site", "Leaf", "Clone", "Treatment", "Replicate")], as.factor)

#######################################################
# Logistic growth model

p <- c(y0 = 1, mumax = 0.3, K = 100)

lower <- c(y0 = 0, mumax = 1e-2, K = 30)
upper <- c(y0 = 10, mumax = 0.6, K = 200)

many_log <- all_growthmodels(
  Count1 ~ Day | Site + Leaf + Clone + Replicate,
  data = my_data,
  p = p,
  upper = upper,
  lower = lower,
  FUN = grow_logistic
)

log_results <- as.data.frame(coef(many_log))

log_results <- log_results %>%
  rownames_to_column(var = "id") %>%
  separate(id, into = c("site", "leaf", "clone", "replicate"), sep = ":", convert = TRUE) %>%
  mutate(
    treatment = case_when(
      replicate %in% c(1, 2, 3) ~ "low",
      replicate %in% c(4, 5, 6) ~ "high"
    ),
    logmu = log(mumax + 1),
    sqrtmu = sqrt(mumax),
    negsqrt = -1 / sqrt(mumax)
  )

#######################################################
# K filtering & estimate

k_filtered <- log_results %>%
  filter(K < 199, K > 30.01) %>%
  group_by(treatment) %>%
  summarise(meanK = mean(K),
            seK = sd(K)/sqrt(n())
            )

r_filtered <- log_results %>%
  filter(mumax < 0.599) %>%
  group_by(treatment) %>%
  mutate(treatment = factor(treatment, levels = c("low", "high"))) %>%
  summarise(meanr = mean(mumax),
            ser = sd(mumax)/sqrt(n()))

  #######################################################
  # K bars
    
    # mean K by treatment
    ggplot(k_filtered, aes(x = treatment, y = meanK)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanK-seK, ymax = meanK+seK), width = 0.35) +
      labs(y="K", x ="Temperature") +
      theme_minimal() 
    
    # mean K by site
    K_filt_site <- log_results %>%
      filter(K < 199, K > 30.01) %>%
      group_by(treatment, site) %>%
      summarise(meanK = mean(K),
                seK = sd(K)/sqrt(n()))
    ggplot(K_filt_site, aes(x = site, y = meanK, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanK-seK, ymax = meanK+seK), position = position_dodge(0.5),width = 0.35) +
      labs(y="K", x ="Site") +
      theme_minimal()
    
  #######################################################
  # r bars
    
    # mean r by treatment
    ggplot(r_filtered, aes(x = treatment, y = meanr, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), width = 0.35) +
      labs(y="r", x ="Temperature") +
      scale_fill_manual(values = c("high" = "#e67147", "low" = "#51969c")) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),            
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
        strip.background = element_rect(fill = "gray90", color = NA), 
        plot.caption = element_text(hjust = 0),
        legend.position = "none",
        axis.text = element_text(size = 14)
      )
    
    # mean r by site
    r_filt_site <- log_results %>%
      filter(mumax < 0.599) %>%
      group_by(treatment, site) %>%
      summarise(meanr = mean(mumax),
                ser = sd(mumax)/sqrt(n()))
    ggplot(r_filt_site, aes(x = site, y = meanr, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), position = position_dodge(0.5),width = 0.35) +
      labs(y="r", x ="Site") +
      theme_minimal()
    
        # no temp
        r_filt_site_only <- log_results %>%
          filter(mumax < 0.599) %>%
          group_by(site) %>%
          summarise(meanr = mean(mumax),
                    ser = sd(mumax)/sqrt(n()))
        ggplot(r_filt_site_only, aes(x = site, y = meanr)) +
          geom_bar(stat = "identity", position = "dodge", width = 0.5) +
          geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), position = position_dodge(0.5),width = 0.35) +
          labs(y="r", x ="Site") +
          theme_minimal()
    
    # mean r by leaf
    r_filt_leaf <- log_results %>%
      filter(mumax < 0.599) %>%
      group_by(treatment, site, leaf) %>%
      summarise(meanr = mean(mumax),
                ser = sd(mumax)/sqrt(n())) %>%
      mutate(leaf = as.factor(leaf))
    ggplot(r_filt_leaf, aes(x = leaf, y = meanr, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), position = position_dodge(0.5),width = 0.35) +
      labs(y="r", x ="Leaf") +
      theme_minimal()
    
    # mean r by clone
    r_filt_clone <- log_results %>%
      filter(mumax < 0.599) %>%
      group_by(treatment, site, leaf, clone) %>%
      summarise(meanr = mean(mumax),
                ser = sd(mumax)/sqrt(n())) %>%
      mutate(leaf = as.factor(leaf),
             clone = as.factor(clone),
             treatment = factor(treatment, levels = c("low", "high")))
    
    ggplot(r_filt_clone, aes(x = clone, y = meanr, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), position = position_dodge(0.5),width = 0.35) +
      scale_fill_manual(values = c("high" = "#e67147", "low" = "#51969c")) +
      labs(y="r", x ="Clone") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),            
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
        strip.background = element_rect(fill = "gray90", color = NA), 
        plot.caption = element_text(hjust = 0),
        legend.position = "none",
        axis.text.x = element_text(size = 14)
      )
    

#######################################################
# Growth rate model

log_lm <- lm(mumax ~ site/leaf/clone, log_results)
loglog_lm <- lm(logmu ~ site/leaf/clone, log_results)

model1 <- lmer(negsqrt ~ treatment + (1|site/leaf/clone), data = log_results)

qqp(resid(model1), "norm")

boxcox(lm(mumax ~ treatment, data = log_results))

simres <- simulateResiduals(fittedModel = model1, plot = TRUE)
testResiduals(simres)

model2 <- lmer(negsqrt ~ treatment + (1|leaf/clone), data = log_results)

anova(model2, model1, test = "LRT")

model3 <- lmer(negsqrt ~ treatment + (1|clone), data = log_results)

anova(model3, model2, test = "LRT")
#######################################################
# Example growth rate plot

# Choose two specific models to plot (e.g., first two randomly selected)
plot_indices <- c(sample(1:102, 2))

# Set up plotting window: 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))

# Loop through the selected models
for (i in plot_indices) {
  plot(many_log[[i]])
  
  # Extract and format metadata for annotations
  this_r <- round(log_results$mumax[i], 4)
  this_K <- round(log_results$K[i], 4)
  this_temp <- log_results$Temperature[i]
  this_id <- paste(log_results$Culture[i], log_results$Replicate[i], sep = "-")
  
  # Add title and annotation
  title(main = paste(this_id, this_temp, sep = "_"), cex.main = 1)
  mtext(paste("r:", this_r), side = 3, adj = 0.15, line = -1, cex = 0.75)
  mtext(paste("K:", this_K), side = 3, adj = 0.15, line = -2, cex = 0.75)
}

