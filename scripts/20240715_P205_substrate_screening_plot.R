# Load the packages
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)

create_enzyme_plot <- function(file_path, plot_width, plot_height, label_size, output_path) {
  # Read the Excel file into a dataframe, selecting only the first 9 columns
  df <- read_excel(file_path, range = cell_cols(1:9))
  
  # Remove the 'Filename' column
  df <- df %>% select(-Filename)
  
  # Filter the dataframe to keep only rows where Time == 24
  filtered_df <- df %>% filter(Time == 24)
  
  # Order the dataframe by descending value of 'Bio'
  filtered_df <- filtered_df %>% arrange(desc(Bio))
  
  # Select the top 3 compounds with the highest 'Bio'
  top_3_filtered_df <- filtered_df %>% slice(1:3)
  
  # Adjust the factor levels of 'Compound' to reflect the order
  top_3_filtered_df$Compound <- factor(top_3_filtered_df$Compound, levels = top_3_filtered_df$Compound[order(-top_3_filtered_df$Bio)])
  
  # Reshape the dataframe for plotting
  long_df <- top_3_filtered_df %>%
    pivot_longer(cols = c(Bio, inactive_enzyme), names_to = "Metric", values_to = "Value") %>%
    pivot_longer(cols = c(Stdev1, Stdev3), names_to = "SD_Metric", values_to = "SD_Value")
  
  # Match standard deviations to the corresponding metrics
  long_df <- long_df %>%
    mutate(SD_Value = ifelse(Metric == "Bio" & SD_Metric == "Stdev1", SD_Value,
                             ifelse(Metric == "inactive_enzyme" & SD_Metric == "Stdev3", SD_Value, NA))) %>%
    filter(!is.na(SD_Value))
  
  # Rename metrics
  long_df <- long_df %>%
    mutate(Metric = recode(Metric, "Bio" = "Enzyme", "inactive_enzyme" = "Inactive enzyme"))
  
  p <- ggplot(long_df, aes(x = Compound, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_errorbar(aes(ymin = Value - SD_Value, ymax = Value + SD_Value), 
                  position = position_dodge(0.9), width = 0.25) +
    labs(x = NULL,  # Remove x-axis label
         y = "Removal",
         fill = NULL) +  # Remove the legend title
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.075, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_fill_manual(values = c("Enzyme" = "blue", "Inactive enzyme" = "grey"), labels = c("P205", "Inactivated P205")) +
    theme_minimal() +
    theme(axis.text.x = element_text(color = "black", vjust = 1, size = label_size),
          axis.text.y = element_text(color = "black", size = label_size),
          axis.title.y = element_text(color = "black", size = label_size, margin = margin(r = 10)),
          plot.title = element_text(color = "black", size = label_size),
          legend.text = element_text(size = label_size),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(0.65, 0.85), 
          legend.direction = "horizontal", 
          legend.background = element_blank(),  # Remove the box around the legend
          legend.box.background = element_blank(),  # Remove the box around the legend
          axis.line = element_line(color = "black"))  # Add axis lines in black
  
  # Save the plot as a PNG file
  ggsave(output_path, plot = p, width = plot_width, height = plot_height, units = "cm")
}

# Example usage
create_enzyme_plot(
  file_path = "data/PLA55254_P205_subtrate_specificity_screening.xlsx",
  plot_width = 20,
  plot_height = 20,
  label_size = 18,
  output_path = "output/PLA_substrate_specificity_plot.png"
)

