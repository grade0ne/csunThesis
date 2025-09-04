data <- read.csv("data/PP_succession_censusdata_cleaned.csv")
data$Bodo <- as.integer(data$Bodo)

protists <- c("Bodo", "Colpoda", 
              "Chrysomonads", "Colpidium", "Cryptomonads", "Cyclidium")


data <- data %>%
  mutate(plantID = as.factor(paste(Plant.ID.Alpha, Plant.ID.Numeric, sep = ""))) %>%
  select(plantID, leafAge, Rotifers, Bodo, Colpoda, 
         Chrysomonads, Colpidium, Cryptomonads, Cyclidium) %>%
  pivot_longer(cols = all_of(protists), names_to = "Protist") %>%
  mutate(
    rotifer_present = case_when(is.na(Rotifers) ~ NA, TRUE ~ Rotifers > 0),
    protist_present = case_when(is.na(value) ~ NA, TRUE ~ value > 0),
    coocurrence = rotifer_present & protist_present) %>%
  group_by(plantID, Protist) %>%
  arrange(leafAge, .by_group = TRUE) %>%
  mutate(
    prev_rotifer = lag(rotifer_present),
    prev_protist = lag(protist_present),
    history = case_when(
      !(coocurrence %in% TRUE) ~ NA_character_,
      is.na(prev_rotifer) | is.na(prev_protist) ~ NA_character_,
      !prev_rotifer & !prev_protist ~ "NPO",
      !prev_rotifer & prev_protist ~ "PF",
      prev_rotifer & !prev_protist ~ "RF",
      prev_rotifer & prev_protist ~ "SCO",
      TRUE ~ NA_character_
    ),
    RPratio = ifelse(!is.na(Rotifers) & !is.na(value) & value > 0, Rotifers / value, NA_real_)
  )

graphdata <- data %>%
  group_by(Protist, history) %>%
  summarize(
    mean = mean(RPratio),
    se = sd(RPratio) / sqrt(length(RPratio)),
    n = length(RPratio)
  )

ggplot(graphdata, aes(x = history, y = mean, fill = history)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(stat = "identity", position = position_dodge(), aes(ymin = mean - se, ymax = mean + se)) +
  labs(x = "History Code", y = "Rotifer:Protist Ratio") +
  facet_wrap( ~ Protist)

# linear model?

data <- data %>%
  mutate(logRatio = log(RPratio))
model1 <- lm(logRatio ~ Protist + history, data = data)

# list: [x] 1) create a co-occurrence check (bool) 
#       [x] 2) look for rotifer presence at previous leaf age time point
#             | Possible scores:
#             | - NPO : No previous occurrences
#             | - RF  : Rotifer first
#             | - PF  : Protist first
#             | - SCO : Sustained co-occurrence
#             |         * look at this separately, does it go up or down? 
#             |           (there are likely not many though)
#       [ ] 3) Compare R:P ratio across different protists; look for co-linearity w/ size