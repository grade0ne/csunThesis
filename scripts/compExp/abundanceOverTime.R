library(tidyverse)
theme_alex<-readRDS('theme_alex.rds')

data<-read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(evolvedTemp,currentTemp,competition,species), as.factor)) %>%
  select(treatID,evolvedTemp,currentTemp,competition,species,day,repNum,count) %>%
  filter(competition == TRUE) %>%
  pivot_wider(names_from = species, values_from = count) %>%
  group_by(evolvedTemp, currentTemp, day, treatID) %>%
  summarise(meanProtist=mean(protist),
            seProtist=sd(protist)/sqrt(length(protist)),
            meanRotifer=mean(rotifer),
            seRotifer=sd(rotifer)/sqrt(length(rotifer))) %>%
  arrange(evolvedTemp, currentTemp, treatID, day) %>%
  mutate(day = day + 0.5,
         fullDay = day %% 1 == 0,
         shapeGroup = interaction(fullDay, evolvedTemp),
         highlight_day = case_when(
           day == 3  ~ "3 d",
           day == 5  ~ "5 d",
           day == 11 ~ "11 d",
           day == 15 ~ "15 d",
           TRUE ~ "Other"))

ggplot(data, aes(x=meanRotifer, y=meanProtist, shape=shapeGroup, 
                 linetype=currentTemp, group=treatID, order=day)) +
  geom_point(stat='identity', size=3) +
  geom_path(linewidth=0.7) +
  scale_shape_manual(values = c("TRUE.25" = 16, "FALSE.25" = 1,
                                "TRUE.30" = 17, "FALSE.30" = 2)) +
  labs(x=expression(paste("Rotifers (0.1 ml"^{-1},")")), 
       y=expression(paste("Protists (0.1 ml"^{-1}, ")")),
       linetype="Assay Temp") +
  scale_x_log10() +
  theme_alex +
  theme(text = element_text(size = 15))


ggplot(data, aes(x = meanRotifer, y = meanProtist,
                 shape = evolvedTemp,          # now shape = evolvedTemp
                 color = highlight_day,        # color highlights by day
                 linetype = currentTemp,
                 group = treatID, order = day)) +
  geom_path(linewidth = 0.7, color = "black") +
  geom_point(stat = "identity", size = 3) +
  scale_shape_manual(
    name = "Evolved Temp",
    values = c("25" = 16,   # filled circle
               "30" = 17)   # filled triangle
  ) +
  scale_color_manual(
    name = "Timepoints",
    values = c(
      "3 d"  = "#1b9e77",
      "5 d"  = "#d95f02",
      "11 d" = "#7570b3",
      "15 d" = "#e7298a",
      "Other"  = "black"
    ),
    breaks = c("3 d", "5 d", "11 d", "15 d")
  ) +
  labs(
    x = expression(paste("Rotifers (0.1 ml"^{-1}, ")")),
    y = expression(paste("Protists (0.1 ml"^{-1}, ")")),
    linetype = "Assay Temp"
  ) +
  #scale_x_log10() +
  theme_alex +
  theme(text = element_text(size = 15))
