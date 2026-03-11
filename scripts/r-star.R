absdata <- read.csv("data/r-star/test1.csv")

library(tidyverse)
library(stringr)

absdata$column <- as.factor(str_remove(absdata$well, "^[A-Za-z]+"))

absdata <- absdata %>%
  mutate(treatment = as.factor(
           case_when(column == 1 ~ "X1",
                               column == 2 ~ "X2",
                               column == 3 ~ "X5",
                               column == 4 ~ "X1mx",
                               column == 6 ~ "medium",
                               column == 7 ~ "DI"
                               )),
         abs = as.numeric(abs)) %>%
  filter(treatment %in% c("X1", "X2", "X5", "X1mx", "medium", "DI"),
         abs < 0.74) %>%
  select(treatment, abs)

abs_summary <- absdata %>%
  group_by(treatment) %>%
  summarize(mean = mean(abs),
            se = sd(abs)/sqrt(length(abs)))
  

ggplot(abs_summary, aes(x = treatment, y = mean)) +
  geom_bar(stat="identity") +
  geom_point(data = absdata, aes(x = treatment, y = abs)) +
  geom_errorbar(stat = "identity", aes(ymin = mean - se, ymax = mean + se))


X1 <- absdata %>%
  filter(treatment == "X1")
hist(X1$abs, breaks = 10)
