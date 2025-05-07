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
    logmu = log(mumax + 1)
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
      summarise(meanK = mean(mumax),
                seK = sd(mumax)/sqrt(n()))
    ggplot(K_filt_site, aes(x = site, y = meanK, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanK-seK, ymax = meanK+seK), position = position_dodge(0.5),width = 0.35) +
      labs(y="K", x ="Site") +
      theme_minimal()
    
  #######################################################
  # r bars
    
    # mean r by treatment
    ggplot(r_filtered, aes(x = treatment, y = meanr)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), width = 0.35) +
      labs(y="r", x ="Temperature") +
      theme_minimal()
    
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
             clone = as.factor(clone))
    ggplot(r_filt_clone, aes(x = clone, y = meanr, fill = treatment)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.5) +
      geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), position = position_dodge(0.5),width = 0.35) +
      labs(y="r", x ="Clone") +
      theme_minimal()

#######################################################
# Growth rate model

log_lm <- lm(mumax ~ site/leaf/clone, log_results)
loglog_lm <- lm(logmu ~ site/leaf/clone, log_results)

model1 <- lmer(mumax ~ treatment + (1|site/leaf/clone), data = log_results)
