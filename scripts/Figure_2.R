#Install packages
pacman::p_load("tidyverse", "ggtree", "DECIPHER", "treeio", "phytools","ape", "stringr", "RColorBrewer", "magrittr")


tr1 <- read.tree("data/trees/Urinary_amidases_controls_unique/20230531_homologs_controls_unique_aligned.fasta.treefile")
tr1$tip.label

midtr <- midpoint.root(tr1)
ggtreem <- ggtree(midtr, color = "gray40", layout = "circular", size = 0.5) # branch.length = "none"

# Metadata
meta_all <- read.csv("data/20230515_metadata_376_amidases_hit.csv") %>% 
  dplyr::mutate(label = word(query_name, sep = "_", 3)) %>% 
  dplyr::mutate(type = "homolog") ## add a column that specifies that these are the homologs


##################
### to check why 2 controls disapper: AYE89271.1_amidase_Comamonas_sp. and WP091572988.1_amidase_Oryzisolibacter_propanilivorax -> they're identical to other two controls

### include labeling of known paracetamol amidase controls
controls <- readAAStringSet("data/7_paracetamol_amidases_controls.txt")

controls_df <- data.frame(
  query_name = names(controls),
  seq = as.character(controls)) %>% 
  dplyr::mutate(label = word(query_name, sep = "_", 1)) %>% # create the label column for merging
  dplyr::mutate(genus = word(query_name, sep = "_", -2)) %>%
  dplyr::mutate(species = word(query_name, sep = "_", -1)) %>%
  mutate(label = ifelse(label == "4YJ6", "4YJ6_1_Chain", label)) %>%
  mutate(label = ifelse(label == "WP091572988.1", "WP_091572988.1", label)) %>% # change name so that it matches the one in meta_all
  dplyr::mutate(type = "paracetamol amidases") # add a column that specifies that these are the controls

### updated hits with correct label
active_homologs <- c("PKZ01497.1", "PKZ94153.1", "PKZ20744.1", "PLB82905.1", "PLB04707.1", "PMC29638.1","PKY89575.1","PKY92403.1","PKY89520.1","PKY70105.1","PKY79758.1","PKZ37694.1","PKY80781.1") 


meta_combined <- controls_df %>%
  select(query_name, seq, label, genus, species, type) %>%
  bind_rows(
    meta_all %>% select(query_name, seq, label, genus, species, type)
  ) %>%
  mutate(
    hits = ifelse(label %in% active_homologs, "yes", "no")
  ) %>%
  distinct(seq, .keep_all = TRUE) %>%   # remove all duplicate hits
  mutate(label = paste0(label,"_", genus, "_", species)) %>% # so that it matches the label of the tree tips
  remove_rownames()

### reading in again the sequences selected for the JGI order
jgi_library <- read.csv("data/selection/JGI_TM1_38_urinary_amidases_order.csv", sep = ";")

### read in CLEAN predictions ###########################################################

### EC NUMBERS OF THE PARACETAMOL AMIDASE HOMOLOGS
clean_homologs <- read.csv("data/CLEAN/20230530_unique_amidases_hits_maxsep_edited.csv", sep = ";", header = TRUE) %>%
  dplyr::mutate(across(c(EC1, EC2, EC3), ~sub("EC:", "", .))) %>%  #remove EC: from the entries
  dplyr::mutate(across(c(EC1, EC2, EC3), ~sub("/.*", "", .))) %>% #remove probabilities including /
  dplyr::mutate(across(c(EC1, EC2, EC3), ~na_if(., ""))) %>% # Add NA to empty cells
  dplyr::rename(label=Identifier) #%>% # rename to match meta dataframe
#dplyr::mutate(label = (word(label, sep = "_", 1))) 

### EC NUMBERS OF THE KNOWN PARACETAMOL AMIDASES (CONTROLS)
clean_controls <- read.csv("data/CLEAN/7_paracetamol_amidases_controls_maxsep.csv", sep = ";", header = TRUE) %>%
  dplyr::mutate(across(c(EC1, EC2), ~sub("EC:", "", .))) %>%  #remove EC: from the entries
  dplyr::mutate(across(c(EC1, EC2), ~sub("/.*", "", .))) %>% #remove probabilities including /
  dplyr::mutate(across(c(EC1, EC2), ~na_if(., ""))) %>% # Add NA to empty cells
  dplyr::rename(label=Identifier) #%>% # rename to match meta dataframe
#dplyr::mutate(label = (word(label, sep = "_", 1)))


#### COMBINING CLEAN PREDICTIONS OF HOMOLOGS WITH CONTROLS
clean_list <- list(clean_homologs, clean_controls)

# List of all potential columns
all_columns <- c( "EC1", "EC2", "EC3")

# Add missing columns to each dataframe
clean_list <- lapply(clean_list, function(df) {
  missing_cols <- setdiff(all_columns, names(df))
  df[, missing_cols] <- NA
  return(df)
})

# Now you can rbind the dataframes
clean_combined <- do.call(rbind, clean_list) 


combined_df <- inner_join(meta_combined, clean_combined, by = "label")
combined_df <- combined_df %>%
  mutate(species = if_else(label == "ANB41810.1_Ochrobactrum_sp.", "sp. TCC-2", species)) %>%
  mutate(species = if_else(label == "ANS81375.1_Ochrobactrum_sp.", "sp. PP-2", species))  %>%
  mutate(EC1 = as.character(substr(EC1, 1, 3))) %>%
  mutate(EC2 = as.character(substr(EC2, 1, 3))) %>%
  mutate(selection = if_else(seq %in% jgi_library$coding_seq, "yes", "no"))

## merge the control dataframe with the homolog metadata
# Read the mapping file
mapping <- read_xlsx("data/20240718_table_SI_p_notation.xlsx")

# Combine the metadata
met <- data.frame(label = ggtreem$data$label[ggtreem$data$isTip]) %>%
  left_join(., combined_df, by = "label") %>%
  mutate(accession_number = sub("^(.+?)_.*", "\\1", label))

# Select only the necessary columns from mapping
mapping_subset <- mapping %>% select(accession_number, p_notation)

# Join met and mapping_subset using accession_number
met <- left_join(met, mapping_subset, by = "accession_number") %>%
  select(-selection, -hits)  %>%
  mutate(p_notation = case_when(
    accession_number == "ANS81375.1" ~ "Mah",
    accession_number == "ANB41810.1" ~ "TccA",
    accession_number == "4YJ6" ~ "4YJ6",
    accession_number == "AFC37599.1" ~ "AmpA",
    TRUE ~ p_notation  # Keep the existing value if no match
  ))



#### PLOTTING THE TREE ##############################  

ggtreemet <- ggtreem %<+% met # add back metadata

tree <- ggtreemet #%>%
  #collapse(node=271) %>%
  #collapse(node=159)

#tree<-ggtreemet  %>%  collapse(node=271) # removes outgroup including PKZ69904, PKZ19817, PKZ18665 (all Gardnerella vaginalis) and PMC32361 (Lb. gasseri)
#tree<-tree  %>%  collapse(node=159) # removes outgroup containing PLT15682 (Lb. fermentum)

## TREE WITHOUT COLORING

p <- tree +
  geom_tiplab(aes(label = paste0("    ", label)), size = 1.5) +
  geom_tippoint(fill= "white", shape = 21, color = "black", size = 2)  # shape 21 is a filled circle, size can be adjusted
p

#ggsave("20231130_IQTREE_white_tips.png", plot = p, width = 12, height = 7, dpi = 900)

##### TREE WITH COLORED EC CLASSES

color_palette <- c("1.1" = "deepskyblue", "1.2" = "blue" , "1.5" = "black", "3.1" = "darkolivegreen1", "3.4" = "chartreuse" , "3.5" = "chartreuse3", "3.6" = "darkolivegreen", "5.1" = "violet", "6.3" = "red" )

p <- tree +
  #geom_tippoint(aes(subset = selection == "yes"), shape = 21, size = 2, stroke = 1.5, color = "orange") +
  geom_tippoint(aes(fill = factor(EC1)), shape = 21, size = 2) +
  scale_fill_manual(values = color_palette) +
  geom_tiplab(aes(label = paste0("    ", label)), size = 1.5)  +
  labs(fill = "Predicted EC")+
  theme(legend.position = "right",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 10, unit = "pt"))
#+
#geom_nodelab(aes(label=node))

p
#ggsave("20231130_IQTREE_tips_colored_by_EC.png", plot = p, width = 12, height = 7, dpi = 900)


### how many amidases are predicted for the first two EC numbers
counting_EC <- sub("(\\d+\\.\\d+).*", "\\1", clean_combined$EC1)
table(counting_EC)

# Count occurrences of each category in counting_EC
counts <- table(counting_EC)
proportions <- counts / sum(counts)
proportions # 53% are 6.3, 40% are 3.5


### COLOR LABELS OF CONTROLS
p <- tree +
  #geom_tippoint(aes(subset = selection == "yes"), shape = 21, size = 2, stroke = 1.5, color = "orange") +
  geom_tippoint(aes(fill = factor(EC1)), shape = 21, size = 2) +
  scale_fill_manual(values = color_palette) +
  geom_tiplab(aes(label = paste0("    ", label, "    ")), size = 1.5)  +
  geom_tiplab(aes(label = paste0("    ", label, "    "), subset = (type == "paracetamol amidases")), size = 1.5, color = "blue") +  # Highlight specific labels
  labs(fill = "Predicted EC")

p
#ggsave("20231130_IQTREE_tips_colored_by_EC_blue_control_labels.png", plot = p, width = 12, height = 7, dpi = 900)

### ADD ORANGE CIRCLE TO THE ONES WE SELECTED
p <- tree +
  geom_tippoint(aes(subset = selection == "yes"), shape = 21, size = 2, stroke = 1.5, color = "orange") +
  geom_tippoint(aes(fill = factor(EC1)), shape = 21, size = 2) +
  scale_fill_manual(values = color_palette) +
  geom_tiplab(aes(label = paste0("    ", label, "    "), 
                  subset = !is.na(genus) & !is.na(species) & type != "paracetamol amidases"), size = 1.5, color = "black") +
  geom_tiplab(aes(label = paste0("    ", label, "    "), 
                  subset = type == "paracetamol amidases"), size = 1.5, color = "blue") +
  labs(fill = "Predicted EC")

p

#ggsave("new_tips_no_collapsed_nodes_two_EC_blue_controls_selection.png", plot = p, width = 12, height = 7, dpi = 900)


### PLOT ONLY ACCESSION NUMBER
p <- tree +
  geom_tippoint(aes(subset = selection == "yes"), shape = 21, size = 2, stroke = 1.5, color = "orange") +
  geom_tippoint(aes(fill = factor(EC1)), shape = 21, size = 2) +
  scale_fill_manual(values = color_palette) +
  geom_tiplab(aes(label = paste0("    ", sub("^(.+?)_.*", "\\1", label),"    "), 
                  subset = type != "paracetamol amidases"), size = 1.5, color = "black") +
  geom_tiplab(aes(label = paste0("    ", sub("^(.+?)_.*", "\\1", label),"    "), 
                  subset = type == "paracetamol amidases"), size = 1.5, color = "grey30") +
  labs(fill = "Predicted EC") +
  theme(legend.position = "bottom") +  # Move legend to bottom
  guides(fill = guide_legend(nrow = 1))  # Display legend items horizontally

p

#ggsave("tree_accession_label_only.png", plot = p, width = 24, height = 14, dpi = 900)


### COLOR BY SELECTION AND HIT
p <- tree +
  geom_tippoint(aes(subset = selection == "yes"), shape = 21, size = 2, stroke = 1.5, color = "orange") +
  geom_tippoint(aes(fill = hits), shape = 21, size = 2) +
  scale_fill_manual(values = c("yes" = "brown", "no" = "white")) +
  #geom_tiplab(aes(label = paste0("    ", label, "    ")), size = 1.5) +
  theme(legend.position = "none")

p
#ggsave("20240413_paracetamol_amidase_hits_tree_not_tiplabel.png", plot = p, width = 12, height = 7, dpi = 900)


#############
# Assuming met dataframe has been created and merged as previously discussed

p <- tree +
  geom_tippoint(aes(fill = factor(EC1)), shape = 21, size = 2) +
  scale_fill_manual(values = color_palette) +
  geom_tiplab(aes(label = paste0("  ",p_notation,"  "), 
                  subset = type != "paracetamol amidases" & !is.na(p_notation)), size = 2, color = "black") +
  geom_tiplab(aes(label = paste0("  ",p_notation,"  "), 
                  subset = type == "paracetamol amidases" & !is.na(p_notation)), size = 2, color = "purple") +
  labs(fill = "Predicted EC") +
  theme(legend.position = "bottom") +  # Move legend to bottom
  guides(fill = guide_legend(nrow = 1))  # Display legend items horizontally

p
ggsave("output/20240718_Figure_2.png", plot = p, width = 16, height = 9, dpi = 900)
