###########################################
###        Create folders in Taco       ###
###      By Marco A. Aquino-Lopez       ###
###########################################



rm(list=ls())
wkdir="~/GitHub/Bayesian_Carbon_Acc/Data/"
setwd(wkdir)

# read data
Data_full <- read.csv('./Cpeat_final_0415.csv')

# Splitting the dataframe by the 'site' column
site_list <- split(Data_full, Data_full$site)

# Loop through each subset of data
for (site_name in names(site_list)) {
  
  # Extract the desired columns
  subset_data <- site_list[[site_name]][, c("depth", "bulk_density", "carbon", "zeroed_age", "biome")]
  
  # Create directory if it doesn't exist
  dir_name <- paste0("~/Documents/Taco/", site_name)
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
  }
  
  # Construct the file path
  file_path <- paste0(dir_name, "/", site_name, ".csv")
  
  # Save the subset data to the file
  write.csv(subset_data, file_path, row.names = FALSE)
}



