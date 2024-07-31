# Load packages
library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")
library("janitor")
library(ggplot2)
library(viridis)
library(dplyr)
library(zoo)
library(pracma)
library(openxlsx)


enzym <- "paracetamol_library_screen"
folder_path <- file.path("data", enzym)
file_path_template <- list.files(folder_path, pattern = "template", full.names = TRUE)
file_path_template <- file_path_template[!grepl("~", file_path_template)]

# Read in the plate template
temp <- read_excel(file_path_template, skip = 2, range = "B3:M11") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = TRUE) %>%
  t() %>%
  as.vector()

temp <- ifelse(temp == "t7 S146A", "S146A", temp)
temp <- ifelse(temp == "T 7 PLA", "lactob", temp)

# Existing OD400 ranges
replicates_ranges_od400 <- list(
  repA = "B995:CU1919",
  repB = "B29:CU944",
  repC = "B28:CU928"
)


all_data_400 <- list()

# Loop through each replicate and its range
for (rep in names(replicates_ranges_od400)) {
  # Construct file path
  tmafils <- list.files(folder_path, pattern = enzym, full.names = TRUE)
  tmafils <- tmafils[!grepl("~|setup|Bradford|screenshot|split|Template|Tris|endpoint", tmafils)]
  tmarep <- tmafils[grep(rep, tmafils)]
  
  # Read in the data for the current replicate
  tma <- read_excel(tmarep, range = replicates_ranges_od400[[rep]], col_types = c("date", rep("numeric", 97)))
  tma[grepl("OVRFLW", tma)] <- 4.0 # Replace "OVRFLW" with max platereader criteria
  
  # Prepare column names
  newnam <- c("time", "temperature_c", paste0(temp))
  colnames(tma) <- make.unique(newnam)
  
  # Add a column to identify the replicate
  tma$replicate <- rep
  
  # Store the processed data
  all_data_400[[rep]] <- tma
}

# Combine all data into a single dataframe
combined_data_400 <- bind_rows(all_data_400)

##########################################################################################################################

# Function to clean and adjust data
clean_and_adjust_data <- function(combined_data) {
  combined_cleaned <- combined_data %>%
    select(-temperature_c) %>%   # Remove the column for temperature (constant at 37C)
    filter(!grepl("Time", time)) %>% # Removes row if column names are duplicated
    mutate(time = round(as.numeric(difftime(time, as.POSIXct("1899-12-30", tz = "UTC"), units = "mins")), 0)) %>%
    select(-contains("NA.")) %>% # Remove columns generated for NAs during import
    select(time, matches("p|lactob|S146A"), replicate) # remove the standard curve data
  
  # Convert from wide to long format while preserving replicate information
  combined_long <- melt(combined_cleaned, id.vars = c('time', 'replicate')) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(variable = gsub("\\.[[:digit:]]", "", variable))
  
  # Adjust time so that all start from 0 and proceed in multiples of 3 minutes
adjusted_combined_long <- combined_long %>%
  group_by(replicate) %>%
  mutate(
    # Identify and segregate each unique series within the replicate, if necessary
    # For simplicity, this step assumes 'replicate' already uniquely identifies each series
    # Subtract the minimum time from each time value within each series
    time = time - min(time),
    # Ensure subsequent readings are in multiples of 3 minutes
    time = round(time / 3) * 3
  ) %>%
  ungroup()
  
  return(adjusted_combined_long)
}

# Apply the function to both OD400
adjusted_combined_long <- clean_and_adjust_data(combined_data_400)


######### AVERAGE TECHNICAL DUPLICATE #################

# Average duplicate measurements within each replicate
averaged_within_replicate <- adjusted_combined_long %>%
  group_by(replicate, time, variable) %>%
  summarise(average_value = mean(value, na.rm = TRUE)) %>%
  ungroup()

######## AVERAGE BIOLOGICAL TRIPLICATES ################
# Calculate the average and standard deviation of the averaged measurements across replicates
average_and_sd_across_replicates <- averaged_within_replicate %>%
  group_by(time, variable) %>%
  summarise(
    overall_average = mean(average_value, na.rm = TRUE),
    sd_across_replicates = sd(average_value, na.rm = TRUE)
  ) %>%
  ungroup()

dat <- average_and_sd_across_replicates
colnames(dat) <- c("time", "variable", "mean", "sd")

####################################
### CALCULATE MEAN OF INACTIVE P205-S146A CONTROL TO NORMALIZE THE DATA TO IT ###
mean_S146A_data <- dat %>%
  filter(variable == "S146A") %>%
  group_by(time) %>%
  summarise(mean_S146A = mean(mean, na.rm = TRUE),.groups = 'drop')


#### calculate rolling mean over 30 min to smooth out the noise ###
dat_adjusted_S146A <- dat %>%
  left_join(mean_S146A_data, by = "time") %>%
  mutate(
    mean_normalized = mean - mean_S146A,
  ) %>%
  arrange(time) %>%
  group_by(variable) %>%
  mutate(
    rolling_mean = rollapply(mean_normalized, width = 10, FUN = mean, fill = NA, align = "center"), # width 10 takes 10 timepoint so it's over 30min
    rolling_sd = rollapply(sd, width = 10, FUN = mean, fill = NA, align = "center")
  ) %>%
  ungroup()

# Ensure that mean_normalized and sd_normalized for S146A are zero
dat_adjusted_S146A <- dat_adjusted_S146A %>%
  mutate(
    rolling_mean = if_else(variable == "S146A", 0, rolling_mean),
    rolling_sd = if_else(variable == "S146A", 0, rolling_sd)
  )

# Plot the data
gg <- ggplot(dat_adjusted_S146A, aes(x = time, y = rolling_mean, group = variable)) +
  geom_line(aes(color = variable)) + # Line for the mean
  geom_line(data = dat_adjusted_S146A %>% filter(variable == "S146A"), aes(x = time, y = rolling_mean, group = variable), color = "red", size = 1.5) +  # Red line for S146A
  scale_color_viridis(discrete = TRUE, option = "D") + # Discrete color scale for lines
  scale_fill_viridis(discrete = TRUE, option = "D") + # Discrete color scale for shading
  labs(x = "Time (minutes)", y = "Mean ± SD", title = "Substrate screening with paracetamol normalized by S146A hydrolysis (red line)") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right") +
  xlim(0, 2600)

print(gg)

#################################################
###### evaluate the normalize OD at time 2600 ####

full_data <- dat_adjusted_S146A %>%
  filter(time == 2601) %>%
  mutate(variable = recode(variable, "lactob" = "P205"),
         variable = toupper(variable))

# Filter data for S146A
dat_S146A <- full_data %>%
  filter(variable == "S146A")

# Calculate the mean and 2 * sd for S146A
mean_S146A <- dat_S146A$rolling_mean
sd_S146A <- dat_S146A$sd
upper_bound_S146A <- mean_S146A + (2 * sd_S146A)
lower_bound_S146A <- mean_S146A - (2 * sd_S146A)

# Ensure dat_time also has the renaming and capitalization applied
dat_time <- dat_adjusted_S146A %>%
  filter(time == 2601) %>%
  mutate(variable = recode(variable, "lactob" = "P205"),
         variable = toupper(variable))

# Create the plot
gg <- ggplot(full_data, aes(x = variable, y = rolling_mean)) +
  geom_point(size = 3, color = "black") + # Points for the mean
  geom_errorbar(aes(ymin = rolling_mean - (2 * rolling_sd), ymax = rolling_mean + (2 * rolling_sd)), width = 0.2, color = "black") + # Error bars for ±2 rolling_sd
  geom_errorbar(data = dat_S146A, aes(ymin = rolling_mean - (2 * sd_S146A), ymax = rolling_mean + (2 * sd_S146A)), width = 0.2, color = "black") + # Error bars for ±2 sd for S146A
  geom_hline(yintercept = upper_bound_S146A, linetype = "dashed", color = "red", size = 1) + # Upper bound for ±2 sd for S146A
  geom_hline(yintercept = lower_bound_S146A, linetype = "dashed", color = "red", size = 1) + # Lower bound for ±2 sd for S146A
  labs(x = "Enzyme", y = "Normalized OD400") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    axis.text = element_text(color = "black"),  # Set axis text color to black
    axis.title = element_text(color = "black"),  # Set axis title color to black
    strip.text = element_text(size = 10, face = "bold", color = "black"),  # Adjust facet label text size and style
    strip.background = element_rect(fill = "lightblue", color = "grey", size = 0.5),  # Add background to facet labels
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),  # Rotate x-axis labels for better readability
    plot.title = element_text(color = "black", hjust = 0.5)  # Set title color to black and center it
  ) +
  coord_cartesian(ylim = c(NA, 1)) 

# Print the plot
print(gg)

#write.csv(full_data, "data/paracetamol_2600min_data.csv", row.names = FALSE)
#ggsave("output/20240723_paracetamol_plot_supplement.png", plot = gg, width = 10, height = 8, dpi = 300)

### WHICH ARE THE HITS?
# Merge full_data with aggregated_data to get mean_S146A and sd_S146A values
threshold_S146A <- 2 * sd_S146A

# Filter rows where rolling_mean is higher than 2 times the sd of S146A
hits <- full_data %>%
  filter(rolling_mean > threshold_S146A)

# Select relevant columns and rename them
hits <- hits %>%
  select(variable, rolling_mean) %>%
  rename(sample_name = variable, Yield = rolling_mean) %>%
  mutate(peak = "paracetamol", wavelength = 400)

table(hits$sample_name) # 12 hits

write.xlsx(hits, file = "paracetamol_hits.xlsx")

