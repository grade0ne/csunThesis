library(tidyverse)

raw_data <- read.csv("data/compexp/compexp2_data.csv")

library(stringr)

data <- raw_data %>%
  mutate(treatment = str_split_i(tube, "", -1),
         replicate = factor(str_split_i(tube, "", -2)),
         evol_temp = factor(case_when(treatment %in% c("a", "b", "c", "d") ~ "25C", TRUE ~ "30C")),
         expr_temp = factor(case_when(treatment %in% c("a", "c", "e", "g") ~ "25C", TRUE ~ "30C")),
         comp_fact = factor(case_when(treatment %in% c("c", "d", "g", "h") ~ TRUE, TRUE ~ FALSE)),
         spec_pres = factor(case_when(replicate %in% c("9", "0") ~ "protist",
                                      as.integer(replicate) %in% c(1:8) & treatment %in% c("c", "d", "g", "h") ~ "both",
                                      TRUE ~ "rotifer")),
         n_rotifer = rotif_temp,
         n_protist = prot_temp * 6.8) %>%
  select(day, replicate, spec_pres, evol_temp, expr_temp, comp_fact, n_rotifer, n_protist)


# ---------- Plots -------------

# Rotifers
figure2_df <- data %>%
  filter(spec_pres %in% c("rotifer", "both"),
         day == "11") %>%
  group_by(evol_temp, expr_temp, comp_fact) %>%
  mutate(evol_temp = case_when(evol_temp == "25C" ~ "History of 25C",
                               evol_temp == "30C" ~ "History of 30C")) %>%
  summarize(mean = mean(n_rotifer), se = sd(n_rotifer)/sqrt(length(n_rotifer)))

ggplot(figure2_df, aes(x = expr_temp, y = mean, fill = comp_fact)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) + 
  geom_errorbar(stat = "identity", position = position_dodge(width = 0.8), 
                aes(ymin = mean - se, ymax = mean + se), width = 0.4) +
  labs(x = "Experiment Temp", y = expression("Rotifers on day 11 (0.1 "*ml^-1*")"), fill = "Competition") +
  facet_wrap(~evol_temp) +
  scale_color_manual(values = c("#1BB6AF", "#E9A17C")) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )


# Protists
figure3_df <- data %>%
  filter(spec_pres %in% c("protist", "both"),
         day == "11") %>%
  group_by(evol_temp, expr_temp, comp_fact) %>%
  mutate(evol_temp = case_when(evol_temp == "25C" ~ "History of 25C",
                               evol_temp == "30C" ~ "History of 30C")) %>%
  summarize(mean = mean(n_protist), se = sd(n_protist)/sqrt(length(n_protist)))

ggplot(figure3_df, aes(x = expr_temp, y = mean, fill = comp_fact)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) + 
  geom_errorbar(stat = "identity", position = position_dodge(width = 0.8), 
                aes(ymin = mean - se, ymax = mean + se), width = 0.4) +
  labs(x = "Experiment Temp", y = expression("Protists on day 11 (0.1 "*ml^-1*")"), fill = "Competition") +
  facet_wrap(~evol_temp) +
  scale_color_manual(values = c("#1BB6AF", "#E9A17C")) +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank()
  )


# ---------- Linear Models -------------

library(car)
library(moments)

# rotifer d11
lm_df <- data %>%
  filter(spec_pres %in% c("rotifer", "both"),
         day == "11") %>%
  group_by(evol_temp, expr_temp, comp_fact)

model1 <- lm(n_rotifer ~ evol_temp * expr_temp * comp_fact, data = lm_df)

# protist d11
lm_df2 <- data %>%
  filter(spec_pres %in% c("protist", "both"),
         day == "11",
         n_protist > 0) %>%
  group_by(evol_temp, expr_temp, comp_fact)

model2 <- lm(n_protist ~ evol_temp * expr_temp * comp_fact, data = lm_df2)
model3 <- lm(log(n_protist) ~ evol_temp * expr_temp * comp_fact, data = lm_df2)
