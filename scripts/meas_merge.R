library(dplyr)


folder_path <- "C:\\Users\\timet\\Documents\\~terHorst\\imaging\\measurements"
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

merged_data <- lapply(csv_files, function(file) {

  data <- read.csv(file)
  
  # file name
  data$file_name <- basename(file)
  
  # file-specific row numbers
  data$row_number <- seq_len(nrow(data))
  
  return(data)
}) %>% 
  bind_rows()

# preview
head(merged_data)

# export
write.csv(merged_data, file = "merged_data.csv", row.names = FALSE)
