library(tidyverse)
library(growthrates)
library(lme4)
library(lmerTest)



# Data #######################################################

data <- read.csv("data/initialMixedCultures/mixedInitial.csv")

# make sure day is numeric and not a factor!

data[c("culture", "temperature", "replicate")] <- 
  lapply(data[c("culture", "temperature", "replicate")], as.factor)

data <- data %>%
  mutate(
    count1 = count + 1
  )


# Growth curves #############################################

data_rep <- data %>%
  group_by(day, culture, temperature) %>%
  mutate(
    mean = mean(count),
    se = sd(count)/sqrt(n()),
    ci = 1.96 * se
  )

ggplot(data_rep, aes(x = day, y = mean, color = temperature, group = interaction(culture, temperature))) +
  geom_point() +
  geom_line() +
  labs(y = "Rotifer Abundance", x = "Day") +
  scale_x_continuous(breaks = seq(3, 24, by = 6)) +
  scale_color_manual(values = c("25" = "blue2", "30" = "red")) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.2, alpha = 0.6) +
  theme_minimal() +
  facet_wrap(~culture) +
  theme(
    panel.grid = element_blank(),            
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
    strip.background = element_rect(fill = "gray90", color = NA), 
    plot.caption = element_text(hjust = 0)  
  )

# Logistic growth model #############################################

p = c(y0 = 2, mumax = 0.3, K = 100)

lower <- c(y0 = 0, mumax = 1e-2, K = 30)
upper <- c(y0 = 10, mumax = 0.6, K = 200)

many_log <- all_growthmodels(
  count1 ~ day | culture + replicate,
  data = data,
  p = p,
  upper = upper,
  lower = lower,
  FUN = grow_logistic
)

log_results <- as.data.frame(coef(many_log))

log_results <- log_results %>%
  rownames_to_column(var = "id") %>%
    separate(id, into = c("Culture", "Replicate"), sep = ":", convert = TRUE) %>%
  mutate(
    Temperature = case_when(
      Replicate %in% c(1, 2, 3) ~ "25C",
      Replicate %in% c(4, 5, 6) ~ "30C"
    ),
    logr = log(mumax + 1)
  )

r_filtered <- log_results %>%
  filter(mumax < 0.599) %>%
  group_by(Temperature) %>%
  summarise(meanr = mean(mumax),
            ser = sd(mumax)/sqrt(n()))

ggplot(r_filtered, aes(x = Temperature, y = meanr, fill = Temperature)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), width = 0.35) +
  scale_fill_manual(values= c("25C" = "blue4", "30C" = "red3")) +
  labs(y="r", x ="Temperature") +
  theme_minimal()

# Plotting predicted values of logistic growth model ###############

plot_indices <- c(sample(1:30, 15), sample(31:60, 15))

par(mfrow = c(6,5), mar = c(2,2,2,1))

for (i in plot_indices) {
  plot(many_log[[i]])
  this_r <- round(log_results$mumax[i], 4)
  this_K <- round(log_results$K[i], 4)
  this_temp <- log_results$Temperature[i]
  this_id <- paste(log_results$Culture[i], log_results$Replicate[i], sep = "-")
  title(main = paste(this_id, this_temp, sep = "_"), cex.main = 1)
  mtext(paste("r:", this_r), side = 3, adj = 0.15, line = -1, cex = 0.65)
  mtext(paste("K:", this_K), side = 3, adj = 0.15, line = -2, cex = 0.65)
}

# linear mixed model #############################################

model1 <- lmer(mumax ~ Temperature + (1|Culture), data= log_results, REML = FALSE)
summary(model1)
isSingular(model1)

model0 <- lm(mumax ~ Temperature, data = log_results)
summary(model0)
