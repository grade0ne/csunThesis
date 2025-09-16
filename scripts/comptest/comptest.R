packages <- c('tidyverse', 'growthrates')
lapply(packages, library, character.only=TRUE)

data <- read.csv("data/comptest/comptestdata.csv")
data[c('treatment', 'rep', 'sp')] <- lapply(data[c('treatment', 'rep', 'sp')], factor)

data <- data %>% mutate(countmL = count * 6.8)

### All counts bar

allBarDF <- data %>%
  filter(day == 5) %>%
  group_by(treatment, sp) %>%
  summarize(mean = mean(countmL),
            n = length(countmL),
            se = sd(countmL) / sqrt(n))

ggplot(allBarDF, aes(x = treatment, y = mean, fill = sp)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(stat = "identity", position = position_dodge(), 
                aes(ymax = mean + se, ymin = mean - se)) +
  facet_wrap(~sp) +
  theme_minimal()

### Protist bar

pBarDF <- data %>%
  group_by(treatment) %>%
  filter(day == 4, sp == "Protist", treatment != "Ronly" ) %>%
  summarize(mean = mean(countmL),
            se = sd(countmL)/sqrt(length(countmL)))

ggplot(pBarDF, aes(x = treatment, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(stat = "identity", position = position_dodge(),
                aes(ymax = mean + se, ymin = mean - se)) +
  theme_minimal()

### Rotifer bar

rBarDF <- data %>%
  group_by(treatment) %>%
  filter(day == 5, sp == "Rotifer", !treatment %in% c("PonlyDI", "PonlyPL")) %>%
  summarize(mean = mean(countmL),
            se = sd(countmL)/sqrt(length(countmL)))

ggplot(rBarDF, aes(x = treatment, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(stat = "identity", position = position_dodge(),
                aes(ymax = mean + se, ymin = mean - se)) +
  theme_minimal()

### Protist scatter time series

pTime <- data %>%
  group_by(day, treatment, sp) %>%
  filter(sp == "Protist", !treatment %in% c("Ronly")) %>%
  summarize(mean = mean(countmL),
            se = sd(countmL)/sqrt(length(countmL))) %>%
  mutate(source = "series")

dayZero <- read.csv("data/comptest/dayZero.csv") %>%
  filter(!treatment %in% c("Ronly")) %>%
  mutate(treatment = factor(treatment),
         se = NA,
         source = "initial")

bridgeData <- rbind(
  pTime[, c("day", "mean", "treatment", "se", "source")],
  dayZero[, c("day", "mean", "treatment", "se", "source")]
)

link_df <- bind_rows(
  dayZero %>% select(day, mean, treatment),
  pTime %>% filter(day == 1) %>% select(day, mean, treatment)
) %>%
  arrange(treatment, day)


pd <- position_dodge(width = 0.15)
ggplot(bridgeData, aes(x = day, y = mean, color = treatment)) +
  annotate("rect", xmin = -Inf, xmax = 0.95, ymin = -Inf, ymax = Inf, 
    fill = "grey95", alpha = 0.5) +
  geom_point(data = bridgeData, stat = "identity", position = pd,
             aes(x = day, y = mean)) +
  geom_errorbar(stat = "identity", position = pd,
                aes(ymax = mean + se, ymin = mean- se)) +
  geom_line(aes(linetype = source), position = pd) +
  geom_point(data = data %>% filter(sp == "Protist", !treatment %in% "Ronly"), 
             aes(group = treatment, color = treatment, x = day, y = countmL), 
             alpha = 0.2, position = pd, size = 1) +
  geom_line(data = link_df, aes(group = treatment, color = treatment), 
            linetype = "dashed", position = pd) +
  
  scale_linetype_manual(values = c("initial" = "dashed", "series" = "solid")) +
  scale_color_manual(values = c("High" = "red2", "Med" = "orange2", "Low" = "green2",
                                "PonlyDI" = "grey60", "PonlyPL" = "black"),
                     labels = c(
                       "High" = "Low comp (p10:1r)",
                       "Med"  = "Med comp (p5:1r)",
                       "Low"  = "High comp (p2:1r)",
                       "PonlyDI" = "Protist only (DI)",
                       "PonlyPL" = "Protist only (media)"),
                     breaks = c("Low", "Med", "High", "PonlyDI", "PonlyPL")) +
  labs(x = "Day", y = "Count (per 0.1 ml)", color = "Treatment") +
  guides(linetype = "none") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, linewidth = .2),
    panel.grid = element_blank()
  )

### Rotifer scatter time series

rTime <- data %>%
  group_by(day, treatment, sp) %>%
  filter(sp == "Rotifer", !treatment %in% c("PonlyDI", "PonlyPL")) %>%
  summarize(mean = mean(countmL),
            se = sd(countmL)/sqrt(length(countmL)))

pd <- position_dodge(width = 0.15)
ggplot(rTime, aes(x = day, y = mean, color = treatment)) +
  geom_point(data = rTime, stat = "identity", position = pd,
             aes(x = day, y = mean)) +
  geom_errorbar(stat = "identity", position = pd,
                aes(ymax = mean + se, ymin = mean - se)) +
  geom_line(position = pd) +
  geom_point(data = data %>% filter(sp == "Rotifer", !treatment %in% c("PonlyPL", "PonlyDI")), 
             aes(group = treatment, color = treatment, x = day, y = countmL), 
             alpha = 0.2, position = pd, size = 1) +
  scale_color_manual(values = c("High" = "red2", "Med" = "orange2", "Low" = "green2",
                                "Ronly" = "black")) +
  labs(x = "Day", y = "Count (per 0.1 ml)", color = "Treatment") +
  guides(linetype = "none") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, linewidth = .2),
    panel.grid = element_blank()
  )
