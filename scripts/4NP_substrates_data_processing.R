### processing the plate reader data of the amidase substrate screening with colorimetric substrates

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
library(scales)
library(openxlsx)

# Make sure to open the .RProj amidase-evolution

folder_layout <- "plate_reader_substrate_screening/layout"
folder_path <- file.path("data", folder_layout)

file_path_templates <- list.files(folder_path, pattern = "TM", full.names = TRUE)
file_path_templates <- file_path_templates[!grepl("~", file_path_templates)]

# Initialize a list to store results (if needed for further processing)
plate_data <- list()
templates <- list()  # Create a new list to store all templates

for (file_path in file_path_templates) {
  temp <- read_excel(file_path, range = "B1:M9") %>%
    janitor::clean_names() %>%
    as.matrix(bycol = TRUE) %>%
    t() %>%
    as.vector()
  
  templates[[basename(file_path)]] <- temp  # Store each template in the list
  plate_data[[basename(file_path)]] <- temp  # Optionally store processed data
}
#####################################
## IMPORT THE 4NP SUBSTRATE PLATE READER DATA
# Define in which range the plate reader data is contained in the export files
ranges_4NP <- list(
  p109_to_p122 = "B41:CU762",
  p124_to_p151 = "B41:CU762",
  p152_to_p168 = "B41:CU762",
  p180_to_p200 ="B29:CU750"
)

all_data_410 <- list()

folder_data <- "plate_reader_substrate_screening"
folder_path <- file.path("data", folder_data)

########TEMPLATE ASSIGNMENT
for (rep in names(ranges_4NP)) {
  tmafils <- list.files(folder_path, full.names = TRUE)
  tmafils <- tmafils[!grepl("~|layout|flutamide", tmafils)]
  templates <- templates[!grepl("~|flutamide", names(templates))]
  tmarep <- tmafils[str_detect(tmafils, rep)]
  
  if (length(tmarep) > 0) {
    tma <- read_excel(tmarep, range = ranges_4NP[[rep]], col_types = c("date", rep("numeric", 97)))
    
    # Select the appropriate template based on the filename or other logic
    template <- templates[grepl(rep, names(templates))]
    
    
    newnam <- c("time", "temperature_c", paste0(unlist(template)))
    colnames(tma) <- make.unique(newnam)
    
    tma$replicate <- rep
    all_data_410[[rep]] <- tma
  }
}


# Combine all data into a single dataframe
combined_data_410 <- bind_rows(all_data_410)


###############
# Function to clean and adjust data

clean_and_adjust_data <- function(combined_data) {
  combined_cleaned <- combined_data %>%
    select(-temperature_c) %>%   # Remove the column for temperature (constant at 37C)
    filter(!grepl("Time", time)) %>% # Removes row if column names are duplicated
    mutate(time = round(as.numeric(difftime(time, as.POSIXct("1899-12-30", tz = "UTC"), units = "mins")), 0)) %>%
    select(-contains("NA.")) %>% # Remove columns generated for NAs during import
    select(-contains("NA"))
    #select(time, contains("p"), replicate) # Keep only time, columns containing "p", and replicate
  
  # Convert from wide to long format while preserving replicate information
  combined_long <- melt(combined_cleaned, id.vars = c('time', 'replicate')) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(variable = gsub("\\.[[:digit:]]", "", variable)) %>%
    mutate(
      substrate = str_extract(variable, "(?<=_).*$"),  # Extract everything after '_'
      variable = gsub("_.*", "", variable)  # Remove everything after and including '_'
    ) %>%
    mutate(substrate = if_else(is.na(substrate), "4-NP", substrate))
  
  combined_long <- combined_long %>%
    filter(!is.na(value))
  
  # Adjust time so that all start from 0 and proceed in multiples of 1 minutes
  adjusted_combined_long <- combined_long %>%
    group_by(replicate) %>%
    mutate(
      # Subtract the minimum time from each time value within each series
      time = time - min(time)
    ) %>%
    ungroup()
  
  return(adjusted_combined_long)
}
# Apply the function
adjusted_combined_long_410 <- clean_and_adjust_data(combined_data_410) %>%
  rename(set = replicate)

######## AVERAGE BIOLOGICAL TRIPLICATES ################

# Calculate the average and standard deviation of the averaged measurements across replicates
dat <- adjusted_combined_long_410 %>%
  group_by(time, variable, set, substrate) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  ) %>%
  ungroup()

################### PLOTTING #######################
# select 4-NP butyrate and trimethyl to visualize self-hydrolysis as red line OR S146A
relevant_sets <- dat %>%
  filter(variable %in% c("S146A"), substrate %in% c("butyrate", "trimethyl"))

# Generate red shades based on the number of unique sets
num_sets <- length(unique(relevant_sets$set))
red_shades <- scales::brewer_pal(palette = "Reds")(num_sets)

gg <- ggplot(dat, aes(x = time, y = mean, group = variable)) +
  geom_line(aes(color = variable)) +  # Line for the mean of all variables
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = variable), alpha = 0.2) +  # SD shading
  geom_line(data = relevant_sets,
            aes(x = time, y = mean, group = set, color = set), size = 1.5) +  # Red line for 4-NP with different sets
  scale_color_manual(values = c(setNames(red_shades, unique(relevant_sets$set)))) +
  scale_fill_viridis(discrete = TRUE, option = "D") +  # Discrete color scale for shading
  labs(x = "Time (minutes)", y = "OD (410nm)", title = "4-NP substrates screen with 4NP substrate spontaneous hydrolysis as red lines") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right") +
  facet_wrap(~substrate, ncol = 1, scales = "free_y")  # Separate panels for each substrate
gg
ggsave("output/20240427_library_screen_4-NP_substrates_4NP_substrates_red.png", plot = gg, width = 11, height = 8, dpi = 300)


#### PLOT THE CONTROLS #####
dat_filtered <- dat %>%
  filter(substrate %in% c("butyrate", "trimethyl") & variable %in% c("4-NP"))

# Creating the plot with the further filtered data
gg <- ggplot(dat_filtered, aes(x = time, y = mean, group = variable)) +
  geom_line(aes(color = set)) + # Line for the mean
  #geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = variable), alpha = 0.2) + # SD shading
  #geom_line(data = dat_filtered %>% filter(variable == "4-NP"), aes(x = time, y = mean, group = variable), color = "red", size = 1.5) +  # Red line for 4-NP
  scale_color_viridis(discrete = TRUE, option = "D") + # Discrete color scale for lines
  scale_fill_viridis(discrete = TRUE, option = "D") + # Discrete color scale for shading
  labs(x = "Time (minutes)", y = "Mean ± SD", title = "Substrate screening with 4NP substrates - 4-NP self hydrolysis") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right") +
  facet_wrap(~substrate, ncol = 1, scales = "free_y") # Separate panels for butyrate and trimethyl
gg

ggsave("output/20240427_library_screen_4-NP-substrate_self_hydrolysis.png", plot = gg, width = 11, height = 8, dpi = 300)


######## NORMALIZING MEANS #######################

# Extract means for 'S146A' with 'butyrate' and 'trimethyl'
mean_S146A_data <- dat %>%
  filter(variable == "S146A") %>%
  group_by(time, set, substrate) %>%
  summarise(mean_S146A = mean(mean, na.rm = TRUE), .groups = 'drop')

# Join the 4NP means back to the original dataset
dat_adjusted_S146A <- dat %>%
  left_join(mean_S146A_data, by = c("time", "set", "substrate")) %>%
  mutate(mean_normalized = mean - mean_S146A )%>%
filter(substrate %in% c("butyrate", "trimethyl"))

### PLOT NORMALIZED MEANS ############### 

gg<- ggplot(dat_adjusted_S146A, aes(x = time, y = mean_normalized, group = variable)) +
  geom_line(aes(color = variable)) + # Line for the mean
  geom_ribbon(aes(ymin = mean_normalized - sd, ymax = mean_normalized + sd, fill = variable), alpha = 0.2) + # SD shading
  geom_line(data = dat_adjusted_S146A %>% filter(variable == "S146A"), aes(x = time, y = mean_normalized, group = variable), color = "red", size = 1.5) +  # Red line for S146A
  scale_color_viridis(discrete = TRUE, option = "D") + # Discrete color scale for lines
  scale_fill_viridis(discrete = TRUE, option = "D") + # Discrete color scale for shading
  labs(x = "Time (minutes)", y = "Mean ± SD", title = "Substrate screening with 4NP substrates normalized by 4NP substrate self hydrolysis (red line)") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right")+
  facet_wrap(~substrate, ncol = 1, scales = "free_y") # Separate panels the substrates
gg
#ggsave("output/20240427_library_screen_4-NP_substrates_substrate_self-hydrolysis_normalized.png", plot = gg, width = 11, height = 8, dpi = 300)

######


### HISTOGRAM
#### PLOT THE MEAN AND STD OF S146A

# Filter the dataset for the time point of interest and relevant variables
dat_time <- dat_adjusted_S146A %>%
  filter(
    (substrate == "butyrate" & time == 100 & variable == "S146A") |
      (substrate == "trimethyl" & time == 600 & variable == "S146A")
  )


# Prepare the full dataset for histogram plotting
full_data <- dat_adjusted_S146A %>%
  filter(
    (substrate == "butyrate" & time == 100) |
      (substrate == "trimethyl" & time == 600)
  ) %>%
  filter(!variable == "4-NP")


# Create the plot using facet_wrap, sd are for each set
gg <- ggplot(full_data, aes(x = mean_normalized, fill = variable)) +
  geom_histogram(binwidth = 0.002, alpha = 0.6, position = "identity") +
  labs(x = "Mean Normalized Value", y = "Count",
       title = "Histogram of Mean Normalized Values at Time = 100 for buyrate and 600 for trimethyl") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~substrate, scales = "free_y") +
  
  # Add vertical lines for mean values for each variable within their specific substrates
  geom_vline(data = dat_time, aes(xintercept = mean_normalized, color = variable),
             linetype = "solid", size = 1) +
  
  # Add dotted lines for +/- 2 SD
  geom_vline(data = dat_time, aes(xintercept = mean_normalized + (2 * sd), color = variable),
             linetype = "dotted", size = 1) +
  geom_vline(data = dat_time, aes(xintercept = mean_normalized - (2 * sd), color = variable),
             linetype = "dotted", size = 1) +
  
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right") +
  coord_cartesian(ylim = c(NA, 2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))  # Automatically set x-axis breaks

# Print the plot
print(gg)

######## split by set and substrate
gg <- ggplot(full_data, aes(x = mean_normalized, fill = variable)) +
  geom_histogram(binwidth = 0.002, alpha = 0.6, position = "identity") +
  labs(x = "Mean Normalized Value", y = "Count",
       title = "Histogram of Mean Normalized Values at Time = 100 for butyrate and 600 for trimethyl") +
  scale_fill_viridis(discrete = TRUE) +
  facet_grid(set ~ substrate, scales = "free_y") +  # Facet by set and substrate
  
  # Add vertical lines for mean values for each variable within their specific substrates
  geom_vline(data = dat_time, aes(xintercept = mean_normalized, color = variable),
             linetype = "solid", size = 1) +
  
  # Add dotted lines for +/- 2 SD
  geom_vline(data = dat_time, aes(xintercept = mean_normalized + (2 * sd), color = variable),
             linetype = "dotted", size = 1) +
  geom_vline(data = dat_time, aes(xintercept = mean_normalized - (2 * sd), color = variable),
             linetype = "dotted", size = 1) +
  
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 10, face = "bold"),  # Adjust facet label text size and style
    strip.background = element_rect(fill = "lightblue", color = "grey", size = 0.5),  # Add background to facet labels
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    panel.spacing = unit(1, "lines")  # Increase space between panels
  ) +
  coord_cartesian(ylim = c(NA, 2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))  # Automatically set x-axis breaks

# Print the plot
print(gg)
#ggsave("output/20240522_histrogram_4-NP_substrates_set_substrate.png", plot = gg, width = 11, height = 8, dpi = 300)

### WHICH ARE THE HITS?

joined_data <- full_data %>%
  inner_join(dat_time, by = c("set", "substrate"), suffix = c("", "_dat_time"))

table(joined_data$sd_dat_time, joined_data$set)

# Filter rows where mean_normalized is higher than 2 times the sd of dat_time
hits <- joined_data %>%
  mutate(threshold = 2 * sd_dat_time) %>%
  filter(mean_normalized > threshold)

# Select relevant columns
hits <- hits %>%
  select(variable, substrate, mean_normalized, threshold) %>%
  rename(sample_name = variable, peak = substrate, Yield = mean_normalized) %>%
  mutate(wavelength = 400)

# Print the hits dataframe
table(hits$peak)

write.xlsx(hits, file = "4NP_substrates_hits.xlsx")


#####
# Prepare the full dataset for histogram plotting and rename variables
library(ggplot2)
library(dplyr)
library(viridis)

# Prepare the full dataset for histogram plotting and rename variables
full_data_butyrate <- dat_adjusted_S146A %>%
  filter(time == 100) %>%
  mutate(variable = recode(variable, "lactob" = "P205",
                           "Cory" = "P203",
                           "Gord" = "P204"),
         variable = toupper(variable)) %>%
  filter(!variable %in% c("BUTYRATE", "TRIMETHYL", "4-NP") & substrate == "butyrate")

full_data_trimethyl <- dat_adjusted_S146A %>%
  filter(time == 600) %>%
  mutate(variable = recode(variable, "lactob" = "P205",
                           "Cory" = "P203",
                           "Gord" = "P204"),
         variable = toupper(variable)) %>%
  filter(!variable %in% c("BUTYRATE", "TRIMETHYL", "4-NP") & substrate == "trimethyl")


# Ensure dat_time also has the renaming and capitalization applied
dat_time_butyrate <- dat_adjusted_S146A %>%
  filter(time == 100) %>%
  mutate(variable = recode(variable, "lactob" = "P205",
                           "Cory" = "P203",
                           "Gord" = "P204"),
         variable = toupper(variable)) %>%
  filter(variable == "S146A" & substrate == "butyrate")

dat_time_trimethyl <- dat_adjusted_S146A %>%
  filter(time == 600) %>%
  mutate(variable = recode(variable, "lactob" = "P205",
                           "Cory" = "P203",
                           "Gord" = "P204"),
         variable = toupper(variable)) %>%
  filter(variable == "S146A" & substrate == "trimethyl")

# Create the plot for flutamide
gg_butyrate <- ggplot(full_data_butyrate, aes(x = variable, y = mean_normalized, fill = variable)) +
  geom_point(size = 3, color = "black") + # Points for the mean_normalized
  geom_errorbar(aes(ymin = mean_normalized - (2 * sd), ymax = mean_normalized + (2 * sd)), width = 0.2, color = "black") + # Error bars for ±2 sd
  labs(x = "Enzyme", y = "Normalized OD410") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ set, scales = "free_x", nrow = 1) +  # Facet by set, one row
  
  # Add red dotted lines for ±2 sd of S146A for each set
  geom_hline(data = dat_time_butyrate, aes(yintercept = mean_normalized + (2 * sd), group = set), linetype = "dotted", color = "red", size = 1) +
  geom_hline(data = dat_time_butyrate, aes(yintercept = mean_normalized - (2 * sd), group = set), linetype = "dotted", color = "red", size = 1) +
  
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    strip.text = element_blank(), # Adjust facet label text size and style
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
    panel.spacing = unit(1, "lines")  # Increase space between panels
  ) +
  coord_cartesian(ylim = c(NA, 2))

# Print the plot
print(gg_butyrate)
#ggsave("output/4NP-butyrate_plot_100.png", plot = gg_butyrate, width = 11, height = 8, dpi = 300)

# Create the plot for nitroacetanilide
gg_trimethyl <- ggplot(full_data_trimethyl, aes(x = variable, y = mean_normalized, fill = variable)) +
  geom_point(size = 3, color = "black") + # Points for the mean_normalized
  geom_errorbar(aes(ymin = mean_normalized - (2 * sd), ymax = mean_normalized + (2 * sd)), width = 0.2, color = "black") + # Error bars for ±2 sd
  labs(x = "Enzyme", y = "Normalized OD410") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ set, scales = "free_x", nrow = 1) +  # Facet by set, one row
  
  # Add red dotted lines for ±2 sd of S146A for each set
  geom_hline(data = dat_time_trimethyl, aes(yintercept = mean_normalized + (2 * sd), group = set), linetype = "dotted", color = "red", size = 1) +
  geom_hline(data = dat_time_trimethyl, aes(yintercept = mean_normalized - (2 * sd), group = set), linetype = "dotted", color = "red", size = 1) +
  
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    strip.text = element_blank(),  # Adjust facet label text size and style
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
    panel.spacing = unit(1, "lines")  # Increase space between panels
  ) +
  coord_cartesian(ylim = c(NA, 2))

# Print the plot
print(gg_trimethyl)
#ggsave("output/4NP-trimethyl_plot_600.png", plot = gg_trimethyl, width = 11, height = 8, dpi = 300)

#### filter the hits ####
# Combine dat_time for both substrates
dat_time_combined <- bind_rows(dat_time_butyrate, dat_time_trimethyl)

# Merge full_data with dat_time_combined to get mean_S146A and sd_S146A values
joined_data <- full_data %>%
  inner_join(dat_time_combined, by = c("set", "substrate"), suffix = c("", "_dat_time"))

# Filter rows where mean_normalized is higher than 2 times the sd of dat_time
hits <- joined_data %>%
  mutate(threshold = 2 * sd_dat_time) %>%
  filter(mean_normalized > threshold)

# Select relevant columns
hits <- hits %>%
  select(variable, substrate, mean_normalized, threshold) %>%
  rename(sample_name = variable, peak = substrate, Yield = mean_normalized) %>%
  mutate(wavelength = 410)

table(hits$peak)

# Save the hits dataset as a CSV file
write.csv(hits, "hits_1250_4NP_substrates.csv", row.names = FALSE)
