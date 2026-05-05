# prerequisites from mendelson_thesis.Rmd

graphdata <- fig2a_data %>%
  filter(currentTemp == "25")

levels(graphdata$competition) <- c("Competitors Absent", "Competitors Present")

graphdata$emmean_bt <- exp(graphdata$emmean)
graphdata$SE_bt <- graphdata$SE * exp(graphdata$emmean)

ggplot(graphdata, aes(x = evolvedTemp, y = emmean_bt, fill = evolvedTemp)) +
  geom_bar(stat="identity", width = 0.6) +
  geom_errorbar(stat="identity", aes(ymin = emmean_bt - SE_bt, ymax = emmean_bt + SE_bt), width = 0.4) +
  labs(y = "Growth Rate (r)", x = "Evolutionary History") +
  scale_fill_manual(values = c("#61aed1", "#ce7a7e")) +
  facet_wrap(~competition) +
  theme_classic() +
  theme(
    legend.position = "none"
  )
