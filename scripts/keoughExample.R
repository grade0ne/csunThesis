# Worked example of triple nested model structure - Mick Keough
# https://mjkeough.github.io/examples/caballes.nb.html#

library(afex)
library(car)
library(lattice)
library(lme4)
library(lmerTest)
library(nlme)
library(VCA)
library(Rmisc)

my_data <- read.csv("data/tempRes2-30d.csv")

