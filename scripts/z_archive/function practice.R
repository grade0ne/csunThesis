f1 <- function(x, y, z) {
  x + y + z
}

# I want: randomly select 5 isolates from each leaf. 
#
# the result should be a table 

data <- read.csv("data/isol1month.csv")

data_true <- subset(data, `week6` == TRUE)

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