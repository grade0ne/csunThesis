library(dplyr)
library(ggplot2)

# Data prep
plot_data <- exp2_growthcurves %>%
  mutate(
    density_per_ml = count * 10,
    # Create a treatment label combining evolved and current temp
    treatment = case_when(
      species == "rotifer" ~ paste0("E", evolvedTemp, "/C", currentTemp),
      species == "protist" ~ paste0("C", currentTemp)
    )
  ) %>%
  group_by(species, treatment, competition, day) %>%
  summarise(
    mean_density = mean(density_per_ml, na.rm = TRUE),
    se_density = sd(density_per_ml, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # 95% CI using t-distribution
    ci_lower = mean_density - qt(0.975, n - 1) * se_density,
    ci_upper = mean_density + qt(0.975, n - 1) * se_density,
    # Truncate lower CI at small positive value for log scale
    ci_lower = pmax(ci_lower, 1)
  )

# Plot
ggplot(plot_data, aes(x = day, y = mean_density,
                      color = treatment,
                      shape = competition,
                      linetype = competition)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.3, alpha = 0.6) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ species, scales = "free_y") +
  scale_y_log10() +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("No competition", "Competition")) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"),
                        labels = c("No competition", "Competition")) +
  labs(
    x = "Day",
    y = "Density (per mL)",
    color = "Treatment",
    shape = "Competition",
    linetype = "Competition"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "italic")
  )
