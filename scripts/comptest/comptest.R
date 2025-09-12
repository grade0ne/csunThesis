packages <- c('tidyverse', 'growthrates')
lapply(packages, library, character.only=TRUE)

data <- read.csv("data/comptest/comptestdata.csv")
data[c('treatment', 'rep', 'sp')] <- lapply(data[c('treatment', 'rep', 'sp')], factor)

graphdata <- data %>%
  filter(day == 2) %>%
  group_by(treatment, sp) %>%
  summarize(mean = mean(count),
            n = length(count),
            se = sd(count) / sqrt(n))

ggplot(graphdata, aes(x = treatment, y = mean, fill = sp)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(stat = "identity", position = position_dodge(), aes(ymax = mean + se, ymin = mean - se)) +
  facet_wrap(~sp) +
  theme_minimal()