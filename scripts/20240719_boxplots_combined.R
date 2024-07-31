# Load necessary libraries
library(dplyr)
library(ggplot2)

df_225nm <- read.csv("data/HPLC/data_225nm.csv")
df_252nm <- read.csv("data/HPLC/data_252nm.csv")
df_305nm <- read.csv("data/HPLC/data_305nm.csv")

# Add a column to distinguish the datasets
df_225nm <- df_225nm %>% mutate(wavelength = 225)
df_252nm <- df_252nm %>% mutate(wavelength = 252)
df_305nm <- df_305nm %>% mutate(wavelength = 305)

# Combine all datasets
df_combined <- bind_rows(df_225nm, df_252nm, df_305nm)

df_filtered <- df_combined %>%
  filter(!is.na(median_value))

average_iqr_plot_filtered <- df_filtered %>%
  ggplot(aes(x = sample_name, y = median_value, color = as.factor(peak))) + 
  geom_point() +
  geom_errorbar(aes(ymin = median_value - 1.5 * iqr_value, ymax = median_value + 1.5 * iqr_value), width = 0.2) +
  geom_hline(data = . %>% filter(sample_name == "S146A"), 
             aes(yintercept = median_value), linetype = "dotted", color = "red") +
  geom_hline(data = . %>% filter(sample_name == "S146A"), 
             aes(yintercept = median_value + 1.5 * iqr_value), linetype = "dotted", color = "blue") +
  geom_hline(data = . %>% filter(sample_name == "S146A"), 
             aes(yintercept = median_value - 1.5 * iqr_value), linetype = "dotted", color = "blue") +
  ylim(-120, 120) +
  geom_hline(yintercept = 50, color = "darkgreen") +
  facet_wrap(~ peak, ncol = 2, scales = "free_x") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text.y = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet background
    axis.title.x = element_blank(),  # Remove individual x-axis labels
    axis.title.y = element_text(size = 12),  # Adjust y-axis label size
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm")
  ) +
  labs(y = "Relative removal (%)",
    color = "Peak"
  )

# Print the plot
print(average_iqr_plot_filtered)
ggsave("boxplot_filtered.png", plot = average_iqr_plot_filtered, width = 12, height = 18, dpi = 600)
