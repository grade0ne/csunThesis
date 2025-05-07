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

my_data <- read.csv("data/tempRes2-30d.csv")

my_data <- my_data %>%
  mutate(
  Count1 = my_data$Count + 1,
  Replicate = as.integer(substr(my_data$ID, nchar(my_data$ID), nchar(my_data$ID))),
  log_count = log(Count1)
)

my_data[c("Site", "Leaf", "Clone", "Treatment", "Replicate")] <-
  lapply(my_data[c("Site", "Leaf", "Clone", "Treatment", "Replicate")], as.factor)


p <- c(y0 = 1, mumax = 0.3, K = 100)

lower <- c(y0 = 0, mumax = 1e-2, K = 30)
upper <- c(y0 = 2, mumax = 0.6, K = 200)

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
    )
  )

mumax_plot <- log_results %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(mumax),
    se = sd(mumax) / sqrt(n())
  )

ggplot(data = mumax_plot, aes(x = treatment, y = mean)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.35)) +
  theme_minimal()

log_lm <- lm(mumax ~ site/leaf/clone, log_results)

model1 <- lmer(mumax ~ treatment + (1|site/leaf/clone), data = log_results)
