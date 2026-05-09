exp1_growth_summary %>%
  group_by(Site, Leaf, Clone) %>%
  summarise(
    mean_mumax = mean(mumax),
    mean_timing = mean(t_mumax),
    .groups = "drop"
  ) %>%
  summarise(cor = cor(mean_mumax, mean_timing))
