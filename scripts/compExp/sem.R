# rotifer r, protist r, rotifer K, protist K, contemporary temp, evolutionary history, competition, rotifer size?
#
#### Data structuring ####

library(tidyverse)
rawData<-read.csv("data/compexp/compexpData.csv") %>%
  mutate(across(c(evolvedTemp,currentTemp,competition,species,treatID,repNum,repID),
                as.factor))

rData <- read.csv('data/compexp/rotiferParameters.csv')

pData <- read.csv('data/compexp/protistParameters.csv')


# Clean and harmonize both datasets

r_df <- rData %>%
  mutate(
    species = "rotifer",
    compID = case_when(
      competition ~ treatment,  
      TRUE ~ str_c("E", evolvedTemp, ifelse(currentTemp == 25, "NR", "HR"))  
    ),
    compKey = case_when(
      competition ~ treatment,
      TRUE ~ ifelse(currentTemp == 25, "HP", "NP") 
    ),
    r_R = mumax,
    K_R = K
  ) %>%
  select(compKey, evolvedTemp, currentTemp, competition, repNum, r_R, K_R)

p_df <- pData %>%
  mutate(
    species = "protist",
    compKey = treatment,
    r_C = mumax,
    K_C = K
  ) %>%
  select(compKey, evolvedTemp, currentTemp, competition, repNum, r_C, K_C)

# Merge by condition key and replicate
sem_wide <- full_join(r_df, p_df,
                      by = c("compKey", "currentTemp", "competition", "repNum"))

# Clean up and add derived variables
sem_wide <- sem_wide %>%
  mutate(
    evolvedTemp = coalesce(evolvedTemp.x, evolvedTemp.y),
    evolvedTemp = as.factor(evolvedTemp),
    currentTemp = as.factor(currentTemp),
    competition = as.logical(competition),
    repNum = as.integer(repNum),
    lnr_R = log(r_R),
    lnr_C = log(r_C),
    competition = factor(competition, levels = c(FALSE, TRUE), labels = c("noComp","withComp")),
    evo_c   = ifelse(evolvedTemp == 30,  +0.5, -0.5),
    assay_c = ifelse(currentTemp == 30,  +0.5, -0.5),
    comp01  = ifelse(competition == "withComp", 1, 0),
    evoXassay = evo_c * assay_c
  ) %>%
  select(evolvedTemp, currentTemp, competition, repNum,
         lnr_R, K_R, lnr_C, K_C, evo_c, assay_c, comp01, evoXassay)

#### SEM ####

library(lavaan)

model_single <- '
  # rotifer performance:
  lnr_R ~ a1*evo_c + a2*assay_c + a3*evoXassay + a4*comp01
  K_R ~ b1*evo_c + b2*assay_c + b3*evoXassay + b4*comp01
  
  # protist performance:
  lnr_C ~ c1*assay_c + c2*comp01
  K_C ~ d1*assay_c + d2*comp01
  
  #covariance
  lnr_R ~~ K_R
  lnr_C ~~ K_C
'

fit_single <- sem(
  model_single,
  data = sem_wide,
  estimator = "MLR",       
  missing   = "fiml",      
  meanstructure = TRUE,
  std.lv = TRUE,
  std.ov = TRUE
)

summary(fit_single, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

library(lavaanPlot)

lbl <- c(
  evo_c   = "Evolved temp",
  assay_c = "Contemporary temp",
  comp01  = "Competition",
  lnr_R   = "Rotifer growth rate (ln r)",
  K_R     = "Rotifer carrying capacity (K)",
  lnr_C   = "Protist growth rate (ln r)",
  K_C     = "Protist carrying capacity (K)",
  evoXassay = "Evolution X Contemp. T"
)

lavaanPlot(
  model = fit_single,
  coefs = TRUE,       
  stand = TRUE,
  sig = 0.05,
  stars = "regress",   
  labels = lbl,
  node_options = list(shape = "box", style = "filled", fillcolor = "gray95",
                      fontname = "Helvetica", color = "grey40"),
  edge_options = list(color = "grey35", penwidth = 1.2),
  graph_options = list(rankdir = "LR")  # left-to-right
)

library(semPlot)

color_vec <- c(
  evo_c = "grey90",
  assay_c = "grey90",
  comp01 = "grey90",
  lnr_R = "lightblue",
  K_R = "lightblue",
  lnr_C = "moccasin",
  K_C = "moccasin"
)

semPaths(
  fit_single,
  what = "std",
  whatLabels = "std",
  nodeLabels = lbl,
  layout = "tree",
  rotation = 2,
  edge.label.cex = 0.9,
  label.cex = 0.9,
  color = list(lat = "white", man = color_vec),  # <- use named color vector
  mar = c(5,5,5,5)
)

