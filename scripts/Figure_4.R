library(ggplot2)
library(caret)

# Load the model from the RDS file
xgb_model_filtered <- readRDS("output/20240709_xgboost_retrain_with_70%_common_features_zerovar_241_removed_features.rds")

# Generate variable importance plot
xgb_imp_filtered <- varImp(xgb_model_filtered, scale = FALSE)
xgb_imp_filtered <- as.data.frame(xgb_imp_filtered$importance)

# Rename the features
rename_features <- function(df) {
  new_names <- c(
    "ANILIDE" = "Aryl amide",
    "MW" = "Molecule molecular mass",
    "MW_bond" = "Amide carbonyl substituent molecular mass",
    "MW_ring" = "Aryl moiety molecular mass",
    "Ncharges" = "Number of charges",
    "N" = "Nitrogen",
    "C" = "Carbon",
    "O" = "Oxygen",
    "R3N" = "Tertiary amine",
    "Cl" = "Chlorine",
    "S" = "Sulfur",
    "R2NH" = "Secondary amine",
    "RCOOR" = "Ester",
    "RINGS" = "Aromatic rings",
    "aa300_size" = "Amino acid 300 size",
    "aa226_size" = "Amino acid 226 size",
    "aa366_size" = "Amino acid 366 size",
    "aa300_secondarystruct" = "Amino acid 300 secondary structure",
    "aa266_size" = "Amino acid 266 size",
    "aa226_polarity" = "Amino acid 266 polarity",
    "aa269_size" = "Amino acid 269 size",
    "aa368_polarity" = "Amino acid 368 polarity",
    "aa329_polarity" = "Amino acid 329 polarity",
    "aa256_polarity" = "Amino acid 256 polarity"
    )
  rownames(df) <- sapply(rownames(df), function(x) ifelse(x %in% names(new_names), new_names[x], x))
  return(df)
}

# Apply renaming
xgb_imp_filtered <- rename_features(xgb_imp_filtered)

# Create the ggplot object
importance_plot <- ggplot(xgb_imp_filtered, aes(x = reorder(rownames(xgb_imp_filtered), Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue", color="black") +
  coord_flip() +
  labs(title = "Chemical and amino acid feature importance",
       x = "Features",
       y = "Importance") +
  theme_classic(base_size = 18) + 
  theme(axis.text.x = element_text(vjust = 1),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18))

# Print the plot
print(importance_plot)


# Save the plot as a high-quality, square image
ggsave("output/20240718_variable_importance_xgboost_plot.png", plot = importance_plot, width = 15, height = 10, dpi = 300)

