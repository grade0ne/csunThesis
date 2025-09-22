
data <- read.csv("data/isol7week.csv")

data_true <- subset(data, X7week == TRUE)

data_split <- split(data_true, data_true$sampleID)

set.seed(123)
selected_isolates <- lapply(data_split, function(x) {
  if (nrow(x)>= 5) {
    x[sample(nrow(x), 5), ]
  } else {
    x
  }
})

result <- do.call(rbind, selected_isolates)

result

length(result$sampleID)