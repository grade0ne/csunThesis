genTime <- groupedData %>%
  group_by(.drop = TRUE) %>%
  filter(species == 'rotifer', competition == 'FALSE', day %in% c(3.5, 8.5)) %>%
  mutate(day = case_when(day == 3.5 ~ "T0", TRUE ~ "T5")) %>%
  select(day, treatID, mean) %>%
  pivot_wider(names_from = day, values_from = mean) %>%
  mutate(g = (log10(T5) - log10(T0)) / log10(2),
         Temp = case_when(treatID %in% c("E25HR", "E25HRP", "E30NR", "E30NRP") ~ "25C",
                          TRUE ~ "30C")) %>%
  group_by(Temp) %>%
  summarise(meanG = mean(g),
            se = sd(g)/sqrt(length(g)))
  

# g = log10(Nt) - log10(N0) / log10