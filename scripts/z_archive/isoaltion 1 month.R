library(tidyverse)

data <- read.csv("data/isol1month.csv")

as.factor(data$site)
as.factor(data$leaf)
as.factor(data$X1month)

data$X1month <- factor(data$X1month, c("TRUE", "LOW","FALSE", "CRASH", "DRY", "SP2"))

data_count_C <- data %>%
  filter(site == "C") %>%
  group_by(leaf, X1month) %>%
  summarize(count = n())

ggplot(data=data_count_C, aes(x=X1month, y=count, fill=X1month)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("TRUE" = "green3", 
                               "LOW" = "gold1", 
                               "FALSE" = "coral2", 
                               "CRASH" = "lightgray", 
                               "DRY" = "gray", 
                               "SP2" = "black")) +
  facet_wrap(~leaf) +
  scale_y_continuous(breaks = 1:15) +
  labs(x = " ", y = "Count", fill = " ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
