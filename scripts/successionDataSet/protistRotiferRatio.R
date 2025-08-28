library(tidyverse)

data<-read.csv("Data/PP_succession_censusdata_cleaned.csv")
data<-data%>%
  mutate(plantID = as.factor(paste(Plant.ID.Alpha, Plant.ID.Numeric)),
         logR = log(Rotifers),
         logP = log(Colpidium),
         RtoP = Rotifers/Colpidium,
         PtoR = Colpidium/Rotifers,
         Cooccur = (Rotifers > 0 & Colpidium > 0))%>%
  filter(Cooccur, Colpidium < 100)%>%
  group_by(plantID)%>%
  select(plantID, leafAge, logR, logP, PtoR, RtoP, Cooccur, Rotifers, Colpidium)

plot(PtoR~leafAge, data = data)
plot(RtoP~leafAge, data = data)
plot(Rotifers~Colpidium, data)
plot(Rotifers~Colpidium, data, xlim = c(0, 40), ylim = c(0, 40))

ggplot(data, aes(x=Colpidium, y=Rotifers, color=plantID)) +
  geom_point(stat="identity", size=2.3, position = position_dodge(.3)) +        
  geom_line() +                                                                 # optional, to show mult. measurs.
  xlim(c(0,40)) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

library(car)
model1<-lm(Rotifers~Colpidium, data = data)                                     # type I regression
qqp(resid(model1), "norm")
plot(model1)
summary(model1)

library(lmodel2)                                                                # type II regression
model2<-lmodel2(Rotifers~Colpidium, range.y="relative", range.x="relative", 
                data=data, nperm=99)

# RMA: m = 0.40354873

data <- data %>%
  ungroup()

ggplot(data, aes(x=Colpidium, y=Rotifers)) +
  geom_point(shape=1) +
  guides(fill="none") + 
  geom_smooth(method="lm", formula=y~x)
