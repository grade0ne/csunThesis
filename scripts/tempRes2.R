library(tidyverse)
library(ggplot)

data <- read.csv("data/tempRes2.csv")

graphData <- data %>%
  filter(Day == 15) %>%
  group_by(Clone, Treatment) %>%
  summarize(mean = mean(Count),
            se = sd(Count)/sqrt(length(Count)),
            site = Site,
            leaf= Leaf)

dodgeWidth <- .75

ggplot(data = graphData, aes(x = Clone, y = mean, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(dodgeWidth), width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(dodgeWidth)) +
  theme_minimal() +
  labs(y = "Rotifers (day 15) +/-SE") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA))
