library(tidyverse)

data <- read.csv("data/experiment 2/size_data/example_size_data_format.csv")

data <- data %>%
  separate(ID, into = c("CultureID", "ImageNum", "PostureCode", "temp1", "temp2", "ROINum"), sep = "_", convert = TRUE)

  
