#clean workspace
rm(list = ls())

#load packages
pacman::p_load("tidyverse", "caret", "rsample", "ranger", "Biostrings", "ggpubr",
               "e1071", "dplyr", "seqinr", "readxl", "fastDummies", "DECIPHER", "xgboost","pROC","writexl", "doParallel","openxlsx")

# load and transform data -------------------------------------------------


# Read in the substrate data 
rawdat <- read_csv("data/machine_learning/272_enzyme_substrate_combinations_20240704.csv") %>%
  dplyr::mutate(substrate = tolower(substrate))

hitdat <- rawdat %>%
  dplyr::group_by(sample_name) %>%
  dplyr::summarise(mean(hit_bool))
hitdat 

substrate <- rawdat %>% 
  dplyr::pull(substrate) %>%
  tolower()
substrate

# Now read in all the chemicals
allchems <- read.csv("data/machine_learning/amides_smiles.csv", sep = ";")

unique(substrate)[unique(substrate) %in% allchems$compound] # the substrates that could be realiably quantified

# Read in chemical properties
df2 <- read_xlsx("data/chemical_descriptors/17_chemical_properties_amidase_substrates_20240703.xlsx")

df3 <- df2 %>%
  dplyr::bind_cols(allchems) %>%
  # dplyr::mutate(smiles = allchems$smiles) %>%
  # left_join(., allchems, by = "smiles") %>%
  dplyr::mutate(compound = allchems$compound) %>%
  dplyr::select(-linear_acyl) %>%
  dplyr::mutate(compound = tolower(compound)) %>%
  dplyr::filter(compound %in% substrate) %>%
  dplyr::mutate(smiles = trimws(smiles))
  #filter(!compound %in% c("butyrate", "trimethyl"))

df3$compound # these are the ones which could be parsed

# read in amino acid sequences, subset for those that have a substrate hit
pnot <- read_xlsx("data/AS_enzymes.xlsx") %>%
  dplyr::mutate(p_notation = case_when(grepl("Gord", p_notation) ~ "Gord", 
                                     TRUE ~ p_notation)) %>%
  dplyr::filter(p_notation %in% rawdat$sample_name) %>%
  arrange(p_notation) %>% #order in ascending p_number
  select(p_notation, seq)

pnot_aa <- AAStringSet(pnot$seq) # convert seq to AA format
names(pnot_aa) <- pnot$p_notation

#Align the sequences
#newaln <- AlignSeqs(pnot_aa)
#BrowseSeqs(newaln)
#writeXStringSet(newaln, "data/alignment/16_enzyme_substrate_combinations_v20240704.fasta")

## identify the start of the Cory enzyme
seqs <- read.alignment("data/alignment/16_enzyme_substrate_combinations_v20240704.fasta", format = "fasta")
tri.pos <- words.pos("grlag",seqs$seq[[1]])
seqs$nam

tri.pos # start of the AS region (according to Lee 2015)
end.pos <- 160  # adjust this to be 160 aa = amidase signature region
#Find it in all sequences
#nuc<-lapply(seqs$seq,function(x) { substr(x,tri.pos,tri.pos+end.pos) }) #also need to offset this
nuc<-lapply(seqs$seq,function(x) { substr(x,tri.pos,tri.pos+end.pos) }) #also need to offset this
nucr<-unlist(nuc)
names(nucr)<-seqs$nam 

appendf <- data.frame(nams = names(nucr), seq = nucr) 

# Check the final dataframe
appendf
#NEW PART TO CONVERT TO NUMERICAL FEATURES 
# Convert amino acids 
source("scripts/R/20240603_convert_seq_5aap.r")
fastatyp <- AAStringSet(appendf$seq)
names(fastatyp) <- appendf$nams
dat_list <- strsplit(as.character(fastatyp), split = "")
names(dat_list) <- names(fastatyp)

# NOTE: the next line will take awhile to run
extract_feat_list <- lapply(1:length(dat_list), function(x) { convert_seq_5aap(dat_list[[x]]) })
sapply(extract_feat_list, length)
extract_feat_list
extract_feat_df <- data.frame(matrix(unlist(extract_feat_list), nrow = length(extract_feat_list), byrow=T), stringsAsFactors=FALSE)
colnames(extract_feat_df) <- paste0("aa", sort(rep(tri.pos:(tri.pos+end.pos), times = 5)), "_", c("polarity", "secondarystruct", "size", "codon_diversity", "charge"))
extract_feat_df$nams <- names(dat_list)
colnames(extract_feat_df)

# Merge with the substrate data

learndat <- extract_feat_df %>%
  left_join(appendf, by = "nams") %>%  # Join with appendf to include motif columns
  left_join(rawdat, by = c("nams" = "sample_name")) %>%
  left_join(df3, by = c("substrate" = "compound"))  %>%
  dplyr::mutate(id = paste0(nams, "_", substrate)) %>%
  dplyr::select(-contains("_charge")) %>% # factor V
  dplyr::select(-substrate) 
  #dplyr::select(-contains("aa"))

colnames(learndat)

###
# Remove variables with nonzero variance (optional)
nozdat <- caret::nearZeroVar(learndat, saveMetrics = TRUE, freqCut = 5, uniqueCut = 2) 
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
which_rem
length(which_rem) # remove 241 features


dat <- learndat %>%
 select(-all_of(which_rem)) %>% 
 mutate(bool = as.factor(learndat$hit_bool)) %>%
 mutate(truth = as.factor(ifelse(bool == 1, "yes", "no"))) %>%
 select(-bool, -hit_bool) %>%
 dplyr::filter(!is.na(smiles)) %>%
 #select(-smiles, -MF)
  select(-smiles, -MF, -seq)
colnames(dat)
# Convert NAs to zero
# dat[is.na(dat)] <- 0
unique(word(dat$id, sep = "_", 2))


# training and test set ---------------------------------------------------
# Set random seed 
set.seed(123456)

# Split into test and training data 
# 80% training
# 20% test
dat_split <- rsample::initial_split(dat, strata = "truth", prop = 0.8)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat) 

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("nams", "truth","id")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "truth","id")]
colnames(x_train)

# Dependent variable
y_train <- dat_train$truth
y_test <- dat_test$truth
y_test # check there is a mix of your two variables


# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$id)

df_test <- data.frame(x_test, stringsAsFactors = F, 
                      row.names = dat_test$id)


# Optional tuning of random forest parameters
mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0), round(ncol(df_train)/1.5))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("gini", "extratrees"),
                       min.node.size = 1)

##################################################################################
################ one RANDOM FOREST ###################################################
##################################################################################
# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger", # change this, see caret vignette
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 3, # increase this to 3 when you run the code 
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy - we're up to 80%
rf$finalModel$prediction.error # out-of-bag error

#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp
ggplot(rf_imp, top = 50) + 
  xlab("") +
  theme_classic()

rf_imp
# Testing set
rf_pred <- predict(rf, newdata = df_test)
rf_pred

cm_rf <- confusionMatrix(rf_pred, as.factor(y_test))
cm_rf
cm_rf[["byClass"]]


########################################################
############ one XGBOOST ###################################
########################################################

# Define the grid for xgbDART
xgb_grid <- expand.grid(
  nrounds = c(50, 100),  # fewer boosting iterations
  max_depth = c(3, 6),
  eta = c(0.1, 0.3),  # learning rate
  gamma = c(0, 1),
  colsample_bytree = c(0.7),
  min_child_weight = c(1),
  subsample = c(0.7),
  rate_drop = c(0.1),  # specific to DART
  skip_drop = c(0.1)   # specific to DART
)


# Train the xgbDART model
xgb_model <- train(
  x = df_train,
  y = y_train,
  method = "xgbDART",
  tuneGrid = xgb_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3, 
                           verboseIter = TRUE, classProbs = TRUE, 
                           savePredictions = "final"),
  verbosity= 0,
  nthread = parallel::detectCores(),  # Use all available CPU cores
  tree_method = "hist"  # Efficient tree construction on CPU
)
# Variable importance for xgbDART
xgb_imp <- varImp(xgb_model, scale = FALSE)
ggplot(xgb_imp, top = 50) + 
  xlab("") +
  theme_classic()

# Testing set prediction for xgbDART
xgb_pred <- predict(xgb_model, newdata = df_test)
cm_xgb <- confusionMatrix(xgb_pred, as.factor(y_test))
cm_xgb[["byClass"]]

# ROC curve for xgbDART
xgb_roc <- pROC::roc(response = ifelse(xgb_model$pred$obs == "no", 0, 1),
                     predictor = ifelse(xgb_model$pred$pred == "no", 0, 1),
                     plot = TRUE)
plot(xgb_roc, type = "s", col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6)


##################################################################################
### RUNNING 10 XGBOOSTS TO OBTAIN MOST ROBUST FEATURES #####
###################################################################################

# Register parallel backend
registerDoParallel(cores = parallel::detectCores())

# Define the function to train the model and extract top features
train_and_extract_features_xgb <- function(seed) {
  set.seed(seed)
  
  # Train the xgbDART model
  xgb_model <- train(
    x = df_train,
    y = y_train,
    method = "xgbDART",
    tuneGrid = xgb_grid,
    trControl = trainControl(method = "repeatedcv", number = 10, repeats = 1, 
                             verboseIter = TRUE, classProbs = TRUE, 
                             savePredictions = "final"),
    verbose = TRUE,
    nthread = parallel::detectCores(),
    tree_method = "hist"
  )
  
  # Variable importance for xgbDART
  xgb_imp <- varImp(xgb_model, scale = FALSE)
  
  # Extract top 50 features and their ranks
  importance_df <- data.frame(
    feature = rownames(xgb_imp$importance),
    importance = xgb_imp$importance[, 1],
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number()) %>%
    head(50)
  
  return(importance_df)
}

# Train the model 10 times with different seeds and extract top 50 features
seeds <- 1:10
all_top_features_xgb <- lapply(seeds, train_and_extract_features_xgb)

# Combine the results into a single data frame
combined_df_xgb <- bind_rows(all_top_features_xgb, .id = "iteration")

# Find the features present in at least 70% of the models
common_features_xgb <- combined_df_xgb %>%
  group_by(feature) %>%
  filter(n() >= 0.7 * length(seeds)) %>%
  summarize(avg_rank = mean(rank), count = n())

# Print the common features and their average ranks
print(common_features_xgb)
# write_csv(common_features_xgb, "output/20240709_10_xgboosts_common_features_nozerovar_241_removed.csv")
common_features_xgb <- read_csv("output/20240709_10_xgboosts_common_features_nozerovar_241_removed.csv")
########################################################
### RUN MODEL ONLY WITH ROBUST FEATURES ##########
########################################################

# Filter the training and test datasets to include only the common features
common_feature_names_xgb <- common_features_xgb$feature

df_train_filtered_xgb <- df_train[, common_feature_names_xgb, drop = FALSE]
df_test_filtered_xgb <- df_test[, common_feature_names_xgb, drop = FALSE]

# Retrain the model using only the common features
xgb_model_filtered <- train(
  x = df_train_filtered_xgb,
  y = y_train,
  method = "xgbDART",
  tuneGrid = xgb_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3, 
                           verboseIter = TRUE, classProbs = TRUE, 
                           savePredictions = "final"),
  verbose = TRUE,
  nthread = parallel::detectCores(),
  tree_method = "hist"
)

# Save the trained model to an RDS file
#saveRDS(xgb_model_filtered, "20240709_xgboost_retrain_with_70%_common_features_zerovar_241_removed_features.rds")

# Load the model from the RDS file
xgb_model_filtered <- readRDS("20240709_xgboost_retrain_with_70%_common_features_zerovar_241_removed_features.rds")

# Generate variable importance plot
xgb_imp_filtered <- varImp(xgb_model_filtered, scale = FALSE)
ggplot(xgb_imp_filtered, top = 50) + 
  xlab("") +
  theme_classic()

# Testing set prediction for xgbDART with common features
xgb_pred_filtered <- predict(xgb_model_filtered, newdata = df_test_filtered_xgb)

# Predictions and confusion matrix
xgb_pred_filtered <- predict(xgb_model_filtered, newdata = df_test_filtered_xgb)
cm_xgb_filtered <- confusionMatrix(xgb_pred_filtered, as.factor(y_test))

# Prepare data for export
feature_importance_xgb <- data.frame(Feature = rownames(xgb_imp_filtered$importance), Overall = xgb_imp_filtered$importance$Overall)

# Prepare the confusion matrix by class with metric names
cm_xgb_filtered_byClass <- as.data.frame(cm_xgb_filtered[["byClass"]])
cm_xgb_filtered_byClass$Metric <- rownames(cm_xgb_filtered_byClass)
cm_xgb_filtered_byClass <- cm_xgb_filtered_byClass[, c("Metric", names(cm_xgb_filtered_byClass)[-length(cm_xgb_filtered_byClass)])]

# Save the ROC plot as an image
roc_plot_path <- "roc_plot.png"
xgb_prob <- predict(xgb_model_filtered, newdata = df_test_filtered_xgb, type = "prob")
xgb_roc <- roc(response = ifelse(y_test == "no", 0, 1), predictor = xgb_prob[, "yes"], plot = TRUE, ci = TRUE)
png(roc_plot_path)
p <-plot(xgb_roc, type = "s", col = "#529DCF", xaxs = "i", yaxs = "i", print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.8, xlim = c(0.5,0), ylim = c(0,1.1))
ggsave(plot=p, filename ="output/ROC_curve_with_CI.png", device ="png", dpi=300)
dev.off()

# Create an Excel workbook and add the feature importance and confusion matrix
wb <- createWorkbook()
addWorksheet(wb, "FeatureImportance")
addWorksheet(wb, "ConfusionMatrixByClass")
addWorksheet(wb, "AUC")

writeData(wb, sheet = "FeatureImportance", feature_importance_xgb)
writeData(wb, sheet = "ConfusionMatrixByClass", cm_xgb_filtered_byClass)

# Insert the ROC plot into the third sheet
insertImage(wb, sheet = "AUC", file = roc_plot_path, width = 6, height = 4, startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, file = "output/20240709_retrained_xgboost_with_70%_common_features_rem241features.xlsx", overwrite = TRUE)

###############################################################
#### RUNNING RANDOM FOREST 10 TIMES TO GET ROBUST FEATURES ####
###############################################################


# Register parallel backend
registerDoParallel(cores = parallel::detectCores())

# Define the function to train the model and extract top features
train_and_extract_features_rf <- function(seed) {
  set.seed(seed)
  
  # Train the Random Forest model
  rf_model <- train(
    x = df_train,
    y = y_train,
    method = "ranger",
    tuneGrid = rf_grid,
    trControl = trainControl(method = "repeatedcv", number = 10, repeats = 1, 
                             verboseIter = TRUE, classProbs = TRUE, 
                             savePredictions = "final"),
    verbose = TRUE,
    importance = "permutation",
    num.threads = parallel::detectCores()
  )
  
  # Variable importance for Random Forest
  rf_imp <- varImp(rf_model, scale = FALSE)
  
  # Extract top 50 features and their ranks
  importance_df <- data.frame(
    feature = rownames(rf_imp$importance),
    importance = rf_imp$importance[, 1],
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number()) %>%
    head(50)
  
  return(importance_df)
}

# Train the model 10 times with different seeds and extract top 50 features
seeds <- 1:10
all_top_features_rf <- lapply(seeds, train_and_extract_features_rf)

# Combine the results into a single data frame
combined_df_rf <- bind_rows(all_top_features_rf, .id = "iteration")

# Find the features present in at least 70% of the models
common_features_rf <- combined_df_rf %>%
  group_by(feature) %>%
  filter(n() >= 0.7 * length(seeds)) %>%
  summarize(avg_rank = mean(rank), count = n())

# Print the common features and their average ranks
print(common_features_rf)
write_csv(common_features_rf, "output/20240709_10_random_forests_common_features_nozerovar_241_removed.csv")

########################################################
### RUN MODEL ONLY WITH ROBUST FEATURES ##########
########################################################

# Filter the training and test datasets to include only the common features
common_feature_names_rf <- common_features_rf$feature

df_train_filtered_rf <- df_train[, common_feature_names_rf, drop = FALSE]
df_test_filtered_rf <- df_test[, common_feature_names_rf, drop = FALSE]

# Retrain the model using only the common features
rf_model_filtered <- train(
  x = df_train_filtered_rf,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3, 
                           verboseIter = TRUE, classProbs = TRUE, 
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation",
  num.threads = parallel::detectCores()
)

# Save the trained model to an RDS file
# saveRDS(rf_model_filtered, "output/20240709_random_forest_retrain_with_70%_common_features_zerovar_241_removed_features.rds")

# Load the model from the RDS file
rf_model_filtered <- readRDS("output/20240709_random_forest_retrain_with_70%_common_features_zerovar_241_removed_features.rds")

# Generate variable importance plot
rf_imp_filtered <- varImp(rf_model_filtered, scale = FALSE)
ggplot(rf_imp_filtered, top = 50) + 
  xlab("") +
  theme_classic()

# Testing set prediction for Random Forest with common features
rf_pred_filtered <- predict(rf_model_filtered, newdata = df_test_filtered_rf)

# Predictions and confusion matrix
rf_pred_filtered <- predict(rf_model_filtered, newdata = df_test_filtered_rf)
cm_rf_filtered <- confusionMatrix(rf_pred_filtered, as.factor(y_test))

# Prepare data for export
feature_importance_rf <- data.frame(Feature = rownames(rf_imp_filtered$importance), Overall = rf_imp_filtered$importance$Overall)

# Prepare the confusion matrix by class with metric names
cm_rf_filtered_byClass <- as.data.frame(cm_rf_filtered[["byClass"]])
cm_rf_filtered_byClass$Metric <- rownames(cm_rf_filtered_byClass)
cm_rf_filtered_byClass <- cm_rf_filtered_byClass[, c("Metric", names(cm_rf_filtered_byClass)[-length(cm_rf_filtered_byClass)])]

# Save the ROC plot as an image
roc_plot_path <- "roc_plot.png"
rf_prob <- predict(rf_model_filtered, newdata = df_test_filtered_rf, type = "prob")
rf_roc <- roc(response = ifelse(y_test == "no", 0, 1), predictor = rf_prob[, "yes"], plot = TRUE)
png(roc_plot_path)
plot(rf_roc, type = "s", col = "#529DCF", xaxs = "i", yaxs = "i", print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6, xlim = c(1.1,-0.1), ylim = c(0,1.1))
dev.off()

# Create an Excel workbook and add the feature importance and confusion matrix
wb <- createWorkbook()
addWorksheet(wb, "FeatureImportance")
addWorksheet(wb, "ConfusionMatrixByClass")
addWorksheet(wb, "AUC")

writeData(wb, sheet = "FeatureImportance", feature_importance_rf)
writeData(wb, sheet = "ConfusionMatrixByClass", cm_rf_filtered_byClass)

# Insert the ROC plot into the third sheet
insertImage(wb, sheet = "AUC", file = roc_plot_path, width = 6, height = 4, startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, file = "output/20240709_retrained_random_forest_with_70%_common_features_rem241features.xlsx", overwrite = TRUE)

