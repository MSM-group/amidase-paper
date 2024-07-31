## alluvial plot with new hits
###############################

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(magrittr)
library(reshape2)
library(tidyverse)
library(viridis)
library(readxl)

dat <- read_excel("data/38_JGI_gene_GMGC_distribution_corrected_16_hits.xlsx") %>%
  select(-'original seq') %>%
  mutate(enzyme = paste0(toupper(shorthand), " [", Genus, "]")) %>%
  rowwise() %>%
  mutate(human_total = sum(c_across(contains("human")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(-contains("human"), human_total) %>%
  rowwise() %>%
  mutate(animal_total = sum(c_across(c("dog-gut", "cat-gut", "mouse-gut", "pig-gut", "animal gut")), na.rm = TRUE)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(aquatic_total = sum(c_across(c("marine", "freshwater")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(-contains("human"), -contains("gut"), -contains("marine"), -contains("freshwater"), human_total, animal_total, aquatic_total) %>%
  select(-other)

dat2 <- dat

# Rename columns
names(dat2)[names(dat2) == "human_total"] <- "human"
names(dat2)[names(dat2) == "animal_total"] <- "animal"
names(dat2)[names(dat2) == "aquatic_total"] <- "aquatic"

# Replace NA entries with 0
dat2[is.na(dat2)] <- 0

  
long_df <- melt(dat2)
################################################################
  
ggplot(data = long_df, 
       aes(y = value, axis1 = enzyme, axis2 = variable)) +
  geom_alluvium(aes(fill = enzyme), alpha = 0.8, width = 0.1) +
  geom_stratum(width = 0.1) +
  geom_text(stat = "stratum", aes(label = paste(after_stat(stratum), after_stat(count), sep = ": ")), size = 5) +
  scale_fill_viridis_d(option = "viridis") +
  theme_void() +
  theme(legend.position = "none") 


##############################################################
# Filter the data to show only alluviums with values greater than 3
long_df_filtered <- long_df %>% filter(value > 3)

# Create the alluvial plot with counts displayed only on the strata
p <-ggplot(data = long_df_filtered, 
       aes(y = value, axis1 = enzyme, axis2 = variable)) +
  geom_alluvium(aes(fill = enzyme), alpha = 0.8, width = 0.1) +
  geom_stratum(width = 0.3) +
  geom_text(stat = "stratum", aes(label = paste(after_stat(stratum), after_stat(count), sep = ": ")), size = 3.5) +
  scale_fill_viridis_d(option = "viridis") +
  theme_void() +
  theme(legend.position = "none")
p

#ggsave("alluvial_plot_high_res.png", plot = p, width = 10, height = 8, dpi = 600)
#############################################################
### now with percentages
p <- ggplot(data = long_df_filtered, 
            aes(y = value, axis1 = enzyme, axis2 = variable)) +
  geom_alluvium(aes(fill = enzyme), alpha = 0.8, width = 0.1) +
  geom_stratum(width = 0.4) +
  geom_text(stat = "stratum", aes(label = paste0(after_stat(stratum), ": ", after_stat(count), " (", round(after_stat(prop) * 100, 1), "%)")), size = 4) +
  scale_fill_viridis_d(option = "viridis") +
  theme_void() +
  theme(legend.position = "none")

# Display the plot
p
# Save the plot as a high-resolution PNG
ggsave("alluvial_plot_high_res.png", plot = p, width = 12, height = 12, dpi = 600)

#### only environmental ########################
dat3 <- dat2 %>%
  select(-human, -animal)

# Convert the wide dataframe to long format using melt
long_df <- melt(dat3, id.vars = c("shorthand", "Genus", "enzyme"))

p <- ggplot(data = long_df, 
            aes(y = value, axis1 = enzyme, axis2 = variable)) +
  geom_alluvium(aes(fill = enzyme), alpha = 0.8, width = 0.1) +
  geom_stratum(width = 0.3) +
  geom_text(stat = "stratum", aes(label = paste0(after_stat(stratum), ": ", after_stat(count), " (", round(after_stat(prop) * 100, 1), "%)")), size = 3.5) +
  scale_fill_viridis_d(option = "viridis") +
  theme_void() +
  theme(legend.position = "none")

# Display the plot
p

ggsave("alluvial_plot_no_human_animal.png", plot = p, width = 12, height = 8, dpi = 600)

#########################################
subset_df <- long_df %>%
  filter(grepl("P131|P148|P162", enzyme))

p <- ggplot(data = subset_df, 
            aes(y = value, axis1 = enzyme, axis2 = variable)) +
  geom_alluvium(aes(fill = enzyme), alpha = 0.8, width = 0.1) +
  geom_stratum(width = 0.4) +
  geom_text(stat = "stratum", aes(label = paste0(after_stat(stratum), ": ", after_stat(count), " (", round(after_stat(prop) * 100, 1), "%)")), size = 3.5) +
  scale_fill_viridis_d(option = "viridis") +
  theme_void() +
  theme(legend.position = "none")

# Display the plot
p

ggsave("alluvial_plot_no_human_animal_promiscuous.png", plot = p, width = 12, height = 8, dpi = 600)
