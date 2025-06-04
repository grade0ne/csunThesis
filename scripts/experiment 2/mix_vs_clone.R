library(tidyverse)
library(growthrates)
library(lme4)
library(lmerTest)

# 1 - Data #####################################################################

data <- read.csv("data/initialMixedCultures/clone_vs_mixed_growth.csv")

data[c("Culture", "Temperature")] <-
  lapply(data[c("Culture", "Temperature")], as.factor)

data <- data %>%
  mutate(
    Count1 = Count + 1
  )

# 2 - Logistic growth model ####################################################

p = c(y0 = 2, mumax = 0.3, K = 100)

lower <- c(y0 = 0, mumax = 1e-2, K = 30)
upper <- c(y0 = 10, mumax = 0.6, K = 200)

many_log <- all_growthmodels(
  Count1 ~ Day | Culture + Temperature + Replicate,
  data = data,
  p = p,
  upper = upper,
  lower = lower,
  FUN = grow_logistic
)

log_results <- as.data.frame(coef(many_log))

log_results <- log_results %>%
  rownames_to_column(var = "id") %>%
  separate(id, into = c("Culture", "Temperature", "Replicate"), sep = ":", convert = TRUE) %>%
  filter(mumax < 0.599) %>%
  mutate(
    Type = as.factor(ifelse(Culture == "Mixed", "Mixed", "Clonal"))
  )

# 3 - Plotting predicted values of logistic growth model #######################
#     * chooses 30 random plots to display. each call of plot_indices assignment           generates a crop of random plots.

plot_indices <- c(sample(1:152, 30))

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

# 4 - Linear mixed model #######################################################

model1 <- lmer(mumax ~ Temperature * Type + (1|Culture), data = log_results)
summary(model1)

graph_data <- log_results %>%
  group_by(Temperature, Type) %>%
  summarise(
    mean_r = mean(mumax),
    se_r = sd(mumax)/sqrt(n()),
    ci_r = 1.96 * se_r,
    mean_K = mean(K),
    se_K = sd(K)/sqrt(n()),
    ci_K = 1.96 * se_K
  )

ggplot(graph_data, aes(x = Temperature, y = mean_r, color = Type)) +
  geom_point(stat = "identity") +
  geom_errorbar(stat = "identity", aes(ymin = mean_r - se_r, ymax = mean_r + se_r), width = 0.15) +
  geom_line(aes(x = Temperature, y = mean_r, group = Type)) +
  labs(y = "Maximum growth rate, r") +
  theme_minimal()

ggplot(graph_data, aes(x = Temperature, y = mean_K, color = Type)) +
  geom_point(stat = "identity") +
  geom_errorbar(stat = "identity", aes(ymin = mean_K - se_K, ymax = mean_K + se_K), width = 0.15) +
  geom_line(aes(x = Temperature, y = mean_K, group = Type)) +
  theme_minimal()
