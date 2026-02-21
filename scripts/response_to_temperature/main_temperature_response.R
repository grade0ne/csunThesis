library(tidyverse)

#-------------- Data Handling --------------

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
  #select(
  #  temperature, site, leaf, clone, area_px) %>%
  mutate(across(c(temperature, site, leaf, clone), factor)) %>%
  filter_out(clone=="6")
  

#-------------- Growth Models --------------

library(growthrates)

INIT_PARAMS <- c(y0 = 1, r = 0.3, alpha = 0.005)
LOWER_PARAMS <- c(y0 = 0, r = 0.001, alpha = 0.001)
UPPER_PARAMS <- c(y0 = 5, r = 2, alpha = 1)

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
out <- grow_logistic_alpha(time, parms = list(y0 = 0, r = 1.1, alpha = 0.02))
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
    table$r2 <- rsquared(fit)
    rownames(table) <- NULL
    
    out_list[[i]] <- table
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

#-------------- Model Building -------------- 
library(car)
library(moments)
library(lme4)
library(lmerTest)

summaryData <- results(models) %>%
  mutate(across(c(Treatment, Site, Leaf, Clone), factor))

summaryData$Tr <- ifelse(summaryData$Treatment == "30C", 0.5, -0.5)

model_raw_r <- lmer(r ~ Tr + (1+Tr || Site:Leaf) + (0+Tr || Site:Leaf:Clone), summaryData)

  vars <- c(
    clone30 = 1.071e-08,
    leaf30 = 6.212e-02,
    leaf25 = 1.969e-02,
    resid = 1.528e-01)
  total_var <- sum(vars)
  
  var_df <- data.frame(
    source = c(names(vars), "total"),
    variance = c(vars, total_var),
    percent = c(vars / total_var * 100, 100),
    row.names = NULL)

model_inv_sqrt_r <- lmer((1/sqrt(r)) ~ Tr + (1+Tr || Site:Leaf) + (0+Tr || Site:Leaf:Clone), summaryData)
  
  hist(1/sqrt(summaryData$r))
  abline(v=mean(1/sqrt(summaryData$r)))
  
  skewness(1/sqrt(summaryData$r))
  kurtosis(1/sqrt(summaryData$r))
  
  qqp(resid(model_inv_sqrt_r), 'norm')
  
  plot(resid(model_inv_sqrt_r)~fitted(model_inv_sqrt_r))
  abline(h=0)

library(emmeans)

emm <- summary(emmeans(model_inv_sqrt_r, ~ Tr))
emm$r <- 1 / (emm$emmean^2)
emm$r_lower <- 1 / (emm$upper.CL^2)    # note: intervals flip
emm$r_upper <- 1 / (emm$lower.CL^2)
emm$r_ci <- emm$r_upper - emm$r

r_lm_result <- emm

sizeData$Tr <- ifelse(sizeData$temperature == "30C", 0.5, -0.5)
  
model_full_size <- lmer(area_px ~ Tr * site + (1+Tr || site:leaf) + (1+Tr || site:leaf:clone), sizeData)

  vars <- c(
    clone30 = 4.283e+03,
    leaf30 = 1.506e-04,
    leaf25 = 4.345e+03,
    resid = 1.217e+06)
  total_var <- sum(vars)
  
  var_df <- data.frame(
    source = c(names(vars), "total"),
    variance = c(vars, total_var),
    percent = c(vars / total_var * 100, 100),
    row.names = NULL)





sizeLeafData <- sizeData %>% group_by(temperature, leaf, clone) %>%
  summarize(mean = mean(area_px), n = length(area_px), se = sd(area_px)/sqrt(length(area_px)))

ggplot(sizeLeafData, aes(x = temperature, y = mean, group=clone, color = leaf)) +
  geom_point(stat = 'identity', position = position_dodge(width = .2)) +
  geom_line(position = position_dodge(width = .2)) +
  geom_errorbar(stat = 'identity', position = position_dodge(width = .2),
                aes(ymin = mean - se, ymax = mean + se)) +
  theme_bw()
  

# ------------ Un-used Models ----------------
  # growth rate
  
  #model0 <- lmer(r ~ Treatment + (1|Site/Leaf/Clone), data = summaryData)
  #
  #model1 <- lmer(log(r) ~ Treatment + (1|Site+Leaf+Clone), data = summaryData)
  #model2 <- lmer((1/r) ~ Treatment + (1|Site/Leaf/Clone), data = summaryData)
  #
  #model2a <- lmer((1/r) ~ Treatment + (1|Site) + (1|Leaf) + (1|Clone), summaryData)
  #
  #model3 <- lm((1/r) ~ Treatment * Clone, summaryData)
  #model3a <- lmer((1/r) ~ Treatment + (1|Clone), summaryData)
  #model3b<- lmer((1/r) ~ Treatment + (1|Leaf/Clone), summaryData)
  #model3c<- lmer((1/r) ~ Treatment + (1|Site/Leaf/Clone), summaryData)
  #model3d<- lmer((1/r) ~ Treatment + (1|Site/Leaf), summaryData)
  #model3e<- lmer((1/r) ~ Treatment + (1|Site), summaryData)
  #
  #model4 <- lm(r ~ Treatment * Clone, summaryData)
  #model4a<- lmer(r ~ Treatment + (1|Clone), summaryData)
  #model4b<- lmer(r ~ Treatment + (1|Leaf/Clone), summaryData)
  #model4c<- lmer(r ~ Treatment + (1|Site/Leaf/Clone), summaryData)
  #model4d<- lmer(r ~ Treatment + (1|Site/Leaf), summaryData)
  #model4e<- lmer(r ~ Treatment + (1|Site), summaryData)
  #
  #m_base  <- lmer(log(r) ~ Treatment + Site + (1|Site:Leaf) + (1|Site:Leaf:Clone), data=summaryData)
  #
  #m_rxn   <- lmer(log(r) ~ Treatment + Site + (1|Site:Leaf) + (1 + Treatment|Site:Leaf:Clone), data=summaryData)
  #
  #model5 <- lmer(log(r) ~ Treatment + (1 | Site) + (1 | Site:Leaf) + 
  #                        (0 + Treatment | Site:Leaf:Clone), data = summaryData)
  #
  #
  #summaryData$Tr <- ifelse(summaryData$Treatment == "30C", 0.5, -0.5)
  ## to remove correlations between random effects using '||', factors have to be replaced by
  ## an explicit contrast-coded predictor
  #
  #summaryData$Treatment <- as.factor(summaryData$Treatment)
  #contrasts(summaryData$Treatment) <- c(-1/2, 1/2)
  #
  #
  #model6 <- lmer(log(r) ~ Tr + (1 + Tr | Site:Leaf:Clone), data = summaryData)
  ## failed to converge, neg. eigenvalue(?) - probably b/c 0 variance in clone @ 25C
  #
  #                                    #*#
  #model7 <- lmer(log(r) ~ Tr + (0 + Tr | Site:Leaf:Clone), data = summaryData)
  ## 
  #
  #
  #model8 <- lmer(r ~ Tr + (1|Site) + (1|Site:Leaf) +
  #              (1 + Tr | Site:Leaf:Clone),
  #              data = summaryData)
  #
  #model9a <- lmer(r ~ Tr + (1 + Tr | Site) + (1 + Tr | Site:Leaf) + (1 + Tr | Site:Leaf:Clone), data = summaryData)
  #
  #model9b <- lmer(r ~ Tr + (1 + Tr | Site:Leaf) + (1 + Tr | Site:Leaf:Clone), data = summaryData)
  #
  #model9e <- lmer(r ~ Treatment + (0 + Treatment | Site:Leaf:Clone), data = summaryData)
  #
  #model9f <- lmer(r ~ Treatment + (1 + Treatment | Site:Leaf) + (1 + Treatment | Site:Leaf:Clone), data = summaryData)
  #
  #model10 <- lmer(r ~ Treatment + Site + (1 + Treatment | Site:Leaf) + (1 + Treatment | Site:Leaf:Clone), summaryData)
  #
  #model11 <- lmer(r ~ Tr + Site + (1+Tr || Site/Leaf) + (1+Tr || Site/Leaf/Clone), summaryData)
  #
  #
  #model12 <- lmer(r ~ Tr + (1+Tr || Site:Leaf) + (1+Tr || Site:Leaf:Clone), summaryData)
  #
  #model13 <- lmer(r ~ Tr + Site + (1+Tr || Site:Leaf) + (0+Tr || Site:Leaf:Clone), summaryData)





## Model Testing:

t1a <- model9
t1b <- model9b # Better fit (x2 = 0, p = 1)
anova(model9, model9b)

t2a <- t1b
t2b <- model9a # Better fit (x2 = 2.38, p = 0.304)
anova(t2a, t2b)

t3a <- t2b
t3b <-

# still not normal w/ ln; gamma link=log?

# size
model9 <- lmer(log(area_px) ~ temperature + (1|site+leaf+clone), data = sizeData)








# ------------ Data Vis ----------------
theme_alex<-readRDS('theme_alex.rds')

## growth rate - total variance: 7.768

# Temperature: F(1,95)=11.009, p=0.0013
ggplot(r_lm_result, aes(x = Tr, y = r)) +
  geom_bar(stat="identity") +
  geom_errorbar(stat="identity", aes(ymin=r-r_lower, ymax=r+r_upper))


graphdata1 <- summaryData %>% group_by(Treatment) %>% 
  summarize(mean=mean(r), se=sd(r)/sqrt(length(r)))
ggplot(graphdata1, aes(x = Treatment, y = mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  labs(y="intrinsic growth rate (r)")

graphdata1a <- summaryData %>% group_by(Treatment, Leaf, Clone) %>%
  summarize(mean=mean(r), se=sd(r)/sqrt(length(r)))
ggplot(graphdata1a, aes(x = Treatment, y = mean, color = Leaf, group = Clone)) +
  geom_point(stat='identity', position = position_dodge(width=0.15)) +
  geom_line(stat='identity', position = position_dodge(width=0.15)) +
  labs(x="Temperature", y="Growth rate") +
  theme_bw() +
  theme(panel.grid = element_blank())


# Site: 1.74% of variance
graphdata2 <- summaryData %>% group_by(Treatment, Site) %>% summarize(mean=mean(r), se=sd(r)/sqrt(length(r)))
ggplot(graphdata2, aes(x = Site, y = mean, fill = Treatment)) +
geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="intrinsic growth rate (r)")

# Leaf: 11.18% of variance
graphdata3 <- summaryData %>% group_by(Treatment, Leaf) %>% summarize(mean=mean(r), se=sd(r)/sqrt(length(r)))
ggplot(graphdata3, aes(x = Leaf, y = mean, fill = Treatment)) +
  geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="intrinsic growth rate (r)")

# Clone: 0% of variance?
graphdata4 <- summaryData %>% group_by(Treatment, Clone) %>% summarize(mean=mean(r), se=sd(r)/sqrt(length(r)))
ggplot(graphdata4, aes(x = Clone, y = mean, fill = Treatment)) +
  geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="intrinsic growth rate (r)")


## size

# Temperature:
graphdata5 <- sizeData %>% group_by(temperature) %>% summarize(mean=mean(area_px), se=sd(area_px)/sqrt(length(area_px)))
ggplot(graphdata5, aes(x = temperature, y = mean, fill = temperature)) +
  geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="body size (px^2)")

# Site:
graphdata6 <- sizeData %>% group_by(temperature, site) %>% summarize(mean=mean(area_px), se=sd(area_px)/sqrt(length(area_px)))
ggplot(graphdata6, aes(x = site, y = mean
                       #, fill = temperature
                       )) +
  geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="body size (px^2)")

# Leaf: 0.03% of variance
graphdata7 <- sizeData %>% group_by(temperature, leaf) %>% summarize(mean=mean(area_px), se=sd(area_px)/sqrt(length(area_px)))
ggplot(graphdata7, aes(x = leaf, y = mean, fill = temperature)) +
  geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="body size (px^2)")

# Clone: 0.72% of variance
graphdata8 <- sizeData %>% group_by(temperature, clone) %>% summarize(mean=mean(area_px), se=sd(area_px)/sqrt(length(area_px)))
ggplot(graphdata8, aes(x = clone, y = mean, fill = temperature)) +
  geom_bar(stat='identity', position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 1), width=.8) +
  labs(y="body size (px^2)")
  

sampleSummary <- sizeData %>%
  group_by(clone) %>% summarize(n25 = sum(temperature=="25C"),
                                n30 = sum(temperature=="30C"))
