#Importing dataset
CommunityMatrix <- read.csv("~/Documents/Input/TabS8.csv")
Metadata <- read.csv("~/Documents/Input/Metadata-TCGA-All-18116-Samples.csv")
#loading depended library 
library(dplyr)
#Removing contaminant genera from CommunityMatrix
#List of genera to remove
ReductGenera <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoeflea", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftia", "Duganella", "Herbaspirillum", 
                  "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobacter", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "Aeromicrobium", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", 
                  "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Paenibacillus", "Streptococcus", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Pedobacter", "Wautersiella", "Deinococcus")
#Iterating through ReductedGenera to remove contaminant genera(ReductGenera) if they exist
for (genus in ReductGenera) {
  if (genus %in% colnames(CommunityMatrix)) {
    CommunityMatrix <- CommunityMatrix %>% select(-intersect(ReductGenera, colnames(CommunityMatrix)))
  }
}

#Applying read threshold of 10 and any value <10 to 0, any NA to 0n 
#CommunityMatrix <- CommunityMatrix %>%
#  mutate(across(everything(), ~ ifelse(is.na(.) | . < 10, 0, .)))

# Apply read threshold of 10 to all genera columns
CommunityMatrix <- CommunityMatrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | . < 10, 0, 1)))

# Filtering metadata to contain only primary tumor samples
colnames(Metadata)[1] <- "KnightID"
filteredMetadata <- Metadata %>%
  filter(sample_type %in% c('Primary Tumor', 'Metastatic')) %>%
  select(KnightID, primary_site, sample_type)
#Combining CommunityMatrix with Metadata for analysis.
Trainingdata <- merge(CommunityMatrix, filteredMetadata, by = "KnightID")
write.csv(Trainingdata, file = "CommunityData.csv", row.names = FALSE)

Trainingdata<- Trainingdata %>% select(-KnightID)

#loading ML dependent libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(randomForest)
library(Boruta)
library(splitstackshape)
library(data.table) # Load data.table for enhanced data manipulation
library(tibble)

#Creating RandomForest Model
#setting random seed for reproducibilty
set.seed(18)
#Removing the few rows with missing values
Trainingdata <- na.omit(Trainingdata)

# Convert all but 'primary_site' and 'sample_type' columns to numeric
Trainingdata <- Trainingdata %>%
  mutate(
    primary_site = as.factor(primary_site),
    sample_type = as.factor(sample_type)
  )

# Filtering data to include only primary tumor and metastatic samples
filtered_data <- Trainingdata %>% filter(sample_type %in% c("Primary Tumor", "Metastatic"))

# Selecting unique cancer types from the 'primary_site' column
cancer_types <- unique(filtered_data$primary_site)

# Initializing lists to store Boruta results and selected features for each cancer type
boruta_results_list <- list()
important_features_list <- list()

# Looping through each cancer type
for (cancer in cancer_types) {
  
  # Subset the data for the current cancer type
  cancer_data <- filtered_data %>% filter(primary_site == cancer)
  
  # Stratified splitting into training and test datasets by 'sample_type'
  set.seed(42) # Set seed for reproducibility
  split_data <- stratified(cancer_data, "sample_type", 0.7, bothSets = TRUE)
  train_data <- split_data$SAMP1
  test_data <- split_data$SAMP2
  
  # Convert train_data to data.table
  setDT(train_data)
  
  # Run Boruta using 'sample_type' as the target (dependent variable) on training data
  boruta_result <- Boruta(sample_type ~ . - primary_site, data = train_data, doTrace = 2, maxRuns = 1500)
  
  # Storing Boruta result in the list
  boruta_results_list[[cancer]] <- boruta_result
  
  # Extracting importance data for plotting
  boruta_importance <- attStats(boruta_result)
  
  # Convert rownames to a column manually using base R
  important_features <- data.frame(feature = rownames(boruta_importance), boruta_importance) %>%
    filter(decision != "Rejected") %>%
    select(feature, meanImp, decision)
  
  # Saving the important features dataframe for comparison
  important_features_list[[cancer]] <- important_features
  
  # Step 4: Plotting Boruta results using ggplot2 and saving as PNG
  ggplot(important_features, aes(x = reorder(feature, meanImp), y = meanImp, fill = decision)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Important Features for", cancer), x = "Features", y = "Importance") +
    theme_minimal() +
    scale_fill_manual(values = c("Confirmed" = "green", "Tentative" = "yellow")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save the plot as PNG
  png_filename <- paste0("Boruta_Feature_Importance_", gsub(" ", "_", cancer), ".png")
  ggsave(png_filename, width = 8, height = 6)
  
  # Select the important columns from the training data using data.table syntax
  selected_features <- important_features$feature
  
  # Ensure selected_features are valid columns in train_data
  important_data <- train_data[, ..selected_features]
  
  # Add back the 'sample_type' column
  important_data[, sample_type := train_data$sample_type]
  
  # Write the important data to a CSV file
  write.csv(important_data, paste0("Selected_Features_", gsub(" ", "_", cancer), ".csv"), row.names = FALSE)
  
  # Print Boruta results and selected features
  print(paste("Boruta output for cancer type:", cancer))
  print(boruta_result)
  print(paste("Selected features for cancer type:", cancer))
  print(selected_features)
}

# Analyzing and comparing the selected features across cancer types



######################################
# Load necessary libraries
library(Boruta)
library(dplyr)
library(data.table)

# Step 1: Ensuring that all columns except for 'sample_type' are numeric and 'sample_type' is a factor
Trainingdata <- Trainingdata %>%
  mutate(
    sample_type = as.factor(sample_type)
  )

# Step 2: Filter the data to include only primary tumor and metastatic samples
filtered_data <- Trainingdata %>% filter(sample_type %in% c("Primary Tumor", "Metastatic"))

# Step 3: Remove the 'primary_site' column, since we're only comparing metastatic vs primary tumor
filtered_data <- filtered_data %>% select(-primary_site)

# Step 4: Convert all microbial features to numeric 

# Step 5: Running Boruta using 'sample_type' as the target (dependent variable)
set.seed(42) # Set seed for reproducibility
boruta_result <- Boruta(sample_type ~ ., data = filtered_data, doTrace = 2)

# Step 6: Print the Boruta output
print("Boruta output for Metastatic vs Primary Tumor comparison")
print(boruta_result)

# Step 7: Plot the Boruta results and save the plot as a PNG file
png_filename <- "Boruta_Feature_Importance_Metastatic_vs_Primary_Tumor.png"
png(file = png_filename, width = 800, height = 600)
plot(boruta_result, cex.axis = 0.7, las = 2, main = "Boruta Feature Importance: Metastatic vs Primary Tumor")
dev.off()

# Step 8: Get the important features selected by Boruta
selected_features <- getSelectedAttributes(boruta_result, withTentative = TRUE)

# Step 9: Print selected features
print("Selected features that distinguish Metastatic vs Primary Tumor:")
print(selected_features)

# Step 10: Save the selected important features data frame for future analysis
# Convert to data.table for easier handling
filtered_data_dt <- as.data.table(filtered_data)

# Select the important columns using .. prefix in data.table
important_data <- filtered_data_dt[, ..selected_features]

# Add back the 'sample_type' column
important_data[, sample_type := filtered_data_dt$sample_type]

# Save the selected important features to a CSV file
write.csv(important_data, "Selected_Features_Metastatic_vs_Primary_Tumor.csv", row.names = FALSE)


