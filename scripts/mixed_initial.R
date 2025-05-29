library(tidyverse)
library(growthrates)

data <- read.csv("data/initialMixedCultures/mixedInitial.csv")

# make sure day is numeric and not a factor!

data[c("culture", "temperature", "replicate")] <- 
  lapply(data[c("culture", "temperature", "replicate")], as.factor)

data <- data %>%
  mutate(
    count1 = count + 1
  )

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
      Replicate %in% c(1, 2, 3) ~ "Low",
      Replicate %in% c(4, 5, 6) ~ "High"
    ),
    logr = log(mumax + 1)
  )

r_filtered <- log_results %>%
  filter(mumax < 0.599) %>%
  group_by(Temperature) %>%
  summarise(meanr = mean(mumax),
            ser = sd(mumax)/sqrt(n()))

ggplot(r_filtered, aes(x = Temperature, y = meanr)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  geom_errorbar(aes(ymin = meanr-ser, ymax = meanr+ser), width = 0.35) +
  labs(y="r", x ="Temperature") +
  theme_minimal()
