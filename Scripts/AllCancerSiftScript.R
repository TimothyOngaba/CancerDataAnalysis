#Reading the data into R
CommunityMatrix <- read.csv("~/Documents/Input/TabS8.csv")
Metadata <- read.csv("~/Documents/Input/Metadata-TCGA-All-18116-Samples.csv")
#loading dependent package
library(dplyr)

#Removing contaminant genera from CommunityMatrix
#List of genera to remove
ReductGenera <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoeflea", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftia", "Duganella", "Herbaspirillum",                "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobacter", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "Aeromicrobium", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Paenibacillus",                   "Streptococcus", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Pedobacter", "Wautersiella", "Deinococcus")
#Iterating through ReductedGenera to remove contaminant genera(ReductGenera) if they exist
for (genus in ReductGenera) {
  if (genus %in% colnames(CommunityMatrix)) {
    CommunityMatrix <- CommunityMatrix %>% select(-intersect(ReductGenera, colnames(CommunityMatrix)))
  }
}

#Applying read threshold of 10 and any value <10 to 0, any NA to 0
CommunityMatrix <- CommunityMatrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | . < 10, 0, .)))

# List of cancer types to analyse
CancerTypes <- c("Brain", "Cervix", "Prostate", "Stomach", "Head and Neck", "Colorectal")

# Filtering metadata to contain only primary tumor samples
colnames(Metadata)[1] <- "KnightID"
filteredMetadata <- Metadata %>%
  filter(sample_type == 'Primary Tumor') %>%
  select(KnightID, primary_site)

# Sorting Metadata respective to cancer type in 'primary_site'
filteredMetadata <- filteredMetadata %>%
  mutate(case = case_when(primary_site %in% CancerTypes ~ primary_site, TRUE ~ 'Other'))

# Reducing metadata to retain only KnightID and primary_site
columns_to_retain <- c("KnightID", "primary_site")
ReducedMetadata <- filteredMetadata %>%
  select(all_of(columns_to_retain))

#Merging CommunityMatrix with Metadata by KnightID the metadata with the raw community matrix
combined_data <- merge(CommunityMatrix, ReducedMetadata, by = "KnightID")

# Initialising list to store the results for each cancer type
cancer_results <- list()

# Loop through each cancer type to perform the analysis
for (cancertype in CancerTypes) {
  # Filter the dataset for the current cancer type
  CancerData <- combined_data %>%
    filter(primary_site == cancertype)
  
  # Filtering combined_data for all other cancer types
  AllOtherCancerTypes <- combined_data %>%
    filter(primary_site != cancertype)
  
  # Initialise results data frame for the current cancer type
  results <- data.frame(Genus = character(), cancer_pos = integer(), cancer_neg = integer(), stringsAsFactors = FALSE)
  
  # Looping through each genus to calculate the number of samples >0 and those equal to zero for the current cancer type
  for (genus in colnames(CancerData)[!colnames(CancerData) %in% "primary_site"]) {
    cancer_pos <- sum(CancerData[[genus]] > 0, na.rm = TRUE)  # Number of samples > 0
    cancer_neg <- sum(CancerData[[genus]] == 0, na.rm = TRUE) # Number of samples == 0
    
    # Adding results to the data frame
    results <- rbind(results, data.frame(Genus = genus, cancer_pos = cancer_pos, cancer_neg = cancer_neg))
  }
  
  # Initialising results data frame for all other cancer types
  result2 <- data.frame(Genus = character(), other_pos = integer(), other_neg = integer(), stringsAsFactors = FALSE)
  
  # Loop through each genus to calculate the number of samples >0 and those equal to zero for other cancer types
  for (genus in colnames(AllOtherCancerTypes)[!colnames(AllOtherCancerTypes) %in% "primary_site"]) {
    other_pos <- sum(AllOtherCancerTypes[[genus]] > 0, na.rm = TRUE)  # Number of samples > 0
    other_neg <- sum(AllOtherCancerTypes[[genus]] == 0, na.rm = TRUE) # Number of samples == 0
    
    # Add the results to the data frame
    result2 <- rbind(result2, data.frame(Genus = genus, other_pos = other_pos, other_neg = other_neg))
  }
  
  # Merging results (current cancer type) with result2 (all other cancers) by Genus
  FinalResults <- merge(results, result2, by = "Genus")
  
  # Function to perform Fisher's exact test and return p-value
  get_p_value <- function(cancer_pos, cancer_neg, other_pos, other_neg) {
    # Create the matrix for Fisher's test
    matrix_data <- matrix(c(cancer_pos, cancer_neg, other_pos, other_neg), nrow = 2, byrow = TRUE)
    
    # Check the matrix dimensions (should be 2x2)
    if (all(dim(matrix_data) == c(2, 2))) {
      test_result <- fisher.test(matrix_data)
      return(test_result$p.value)
    } else {
      # Return NA if the matrix is not 2x2
      return(NA)
    }
  }
  
  # Applying function to each row and creating a new column for p-values
  FinalResults <- FinalResults %>%
    rowwise() %>%
    mutate(p_value = get_p_value(cancer_pos, cancer_neg, other_pos, other_neg)) %>%
    ungroup()
  
  # Sorting FinalResults by p-value in ascending order
  sorted_FinalResults <- FinalResults %>%
    arrange(p_value)
  
  # Adjusting p-values by Benjamini-Hochberg method
  sorted_FinalResults <- sorted_FinalResults %>%
    mutate(p_adjusted = p.adjust(p_value, method = 'BH'))
  
  # Adjust p-values using the Benjamini-Hochberg method
  FinalResults <- FinalResults %>%
    mutate(p_adjusted = p.adjust(p_value, method = 'BH'))
  
  # Identify significantly different genera (p_adjusted < 0.05)
  significant_genera <- FinalResults %>%
    filter(p_adjusted < 0.05)
  
  # Count the number of significantly different genera
  num_significant_genera <- nrow(significant_genera)
  
  # Save the significant genera and their counts
  cancer_results[[cancertype]] <- list(
    FinalResults = FinalResults,
    significant_genera = significant_genera,
    num_significant_genera = num_significant_genera
  )
  
  # Save the sorted results to a CSV file
  output_filename <- paste0("FinalResults_", cancertype, ".csv")
  write.csv(FinalResults, output_filename, row.names = FALSE)
  
  # Print the number of significantly different genera
  print(paste("Number of significantly different genera in", cancertype, ":", num_significant_genera))
}

# Compare the number of significantly different genera across cancer types
significant_genera_counts <- data.frame(
  CancerType = names(cancer_results),
  NumSignificantGenera = sapply(cancer_results, function(x) x$num_significant_genera)
)

# Print the comparison table
print(significant_genera_counts)

# Identifying genera that are highly enriched in each cancer type compared to all others
highly_enriched_genera <- list()

for (cancertype in CancerTypes) {
  enriched_genera <- cancer_results[[cancertype]]$FinalResults %>%
    filter((cancer_pos/(cancer_pos+cancer_neg)) > (other_pos/(other_pos+other_neg))) %>%
    arrange(p_value)
    highly_enriched_genera[[cancertype]] <- enriched_genera
}

# Print highly enriched genera for each cancer type
for (cancertype in CancerTypes) {
  print(paste("Highly enriched genera in", cancertype, ":"))
  print(highly_enriched_genera[[cancertype]])
}

#Plotting number of significantly enriched genera per cancer type
#loadind dependent package
library(ggplot2)

# Create the data frame for plotting
significant_genera_counts <- data.frame(
  CancerType = names(highly_enriched_genera),
  NumSignificantGenera = sapply(highly_enriched_genera, nrow)
)

# Plotting the bar graph
ggplot(significant_genera_counts, aes(x = CancerType, y = NumSignificantGenera)) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(aes(label = NumSignificantGenera), vjust = -0.5) + # Add counts on top of bars
  labs(title = "Number of Significantly Enriched Genera by Tumor Type",
       x = "Tumor Type",
       y = "Number of Significantly Enriched Genera (p < 0.1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##########################
#Perforing PCA on cancer types and genera
# Standardising Combined data
numeric_data <- combined_data %>% select_if(is.numeric)
non_constant_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) > 0]
standardised_data <- scale(non_constant_data)

# Perform PCA
pca_result <- prcomp(standardised_data, center = TRUE, scale. = TRUE)

# Calculating the proportion of variance explained by Principal components
percent_variance <- ((pca_result$sdev^2) / sum(pca_result$sdev^2)) * 100

# Plot PCA
library(ggplot2)
pca_df <- data.frame(pca_result$x, KnightID = combined_data$KnightID, CancerType = combined_data$primary_site)
# Identify samples with the lowest and highest variance for PC2
min_pc2_sample <- pca_df[which.min(pca_df$PC2), ]
max_pc2_sample <- pca_df[which.max(pca_df$PC2), ]

# Print the details of these samples
cat("Sample with the lowest variance for PC2:\n")
print(min_pc2_sample)

cat("Sample with the highest variance for PC2:\n")
print(max_pc2_sample)

# Plot PCA with annotations for the samples with the lowest and highest PC2 values
ggplot(pca_df, aes(x = PC1, y = PC2, color = CancerType)) +
  geom_point(size = 2) +
  geom_text(data = min_pc2_sample, aes(label = paste(KnightID, CancerType, sep = "\n")), hjust = 1.5, vjust = 1.5, color = "blue") +
  geom_text(data = max_pc2_sample, aes(label = paste(KnightID, CancerType, sep = "\n")), hjust = 1.5, vjust = 1.5, color = "red") +
  labs(title = "PCA of Cancer Type Genera",
       x = paste0("PC1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance)"),
       y = paste0("PC2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance)")) +
  theme_minimal()
##########################
# Filter the data for Stomach and Lung samples for further PCA as they have maximum PCA variance
stomach_lung_data <- combined_data %>%
  filter(primary_site %in% c("Stomach", "Lung"))

#Selecting all numeric columns
numeric_data <- stomach_lung_data %>% select_if(is.numeric)

# Remove any constant columns (columns with zero variance)
non_constant_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) > 0]
# Standardizing data
standardised_data <- scale(non_constant_data)
# Performing PCA
pca_result <- prcomp(standardised_data, center = TRUE, scale. = TRUE)
# Create a dataframe for PCA plotting
pca_df <- data.frame(pca_result$x, 
                     KnightID = stomach_lung_data$KnightID, 
                     PrimarySite = stomach_lung_data$primary_site)
# Ploting the PCA results
ggplot(pca_df, aes(x = PC1, y = PC2, color = PrimarySite, label = KnightID)) +
  geom_point(size = 2) +
  labs(title = "PCA of stomach and Lung Samples", 
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) +
  theme_minimal() +
  geom_text(aes(label = KnightID), hjust = 0, vjust = 1.5, size = 3)
# Extracting the loadings (contributions of the variables to the PCs)
loadings <- pca_result$rotation
# Extract loadings for PC1 and PC2
pc1_loadings <- loadings[, 1]
pc2_loadings <- loadings[, 2]

# Sorting the loadings by absolute value to see which variables contribute the most
pc1_contributors <- sort(abs(pc1_loadings), decreasing = TRUE)
pc2_contributors <- sort(abs(pc2_loadings), decreasing = TRUE)

# Viewing the top contributors
print("Top contributors to PC1:")
print(pc1_contributors)

print("Top contributors to PC2:")
print(pc2_contributors)




#############################
#Determining which factors/genera contribute to PC1. Rotation contains the loadings/weights of each PC component
# Loadings of each feature (genus) for the principal components
loadings <- pca_result$rotation

# For the top contributing genera to the first principal component
top_contributors_PC1 <- sort(abs(loadings[,1]), decreasing = TRUE)
top_contributors_PC1 <- head(top_contributors_PC1, n = 10)  

# Print the genera contributing most to PC1
print(top_contributors_PC1)

##########################
#Performing PCoA
# Loading required package
install.packages("vegan")
library(vegan)

#Reading in data and filtering for only neumeric data
genus_data <- combined_data[, -which(colnames(combined_data) == "primary_site")]
genus_data <- genus_data %>% select_if(is.numeric)

# Removing rows with all-zero values in `genus_data`
non_empty_data <- genus_data[rowSums(genus_data != 0, na.rm = TRUE) > 0, ]

# Filtering the metadata correspondingly
filtered_metadata <- combined_data[rowSums(genus_data != 0, na.rm = TRUE) > 0, ]

# Calculate Bray-Curtis distances on non-empty data
bray_dist <- vegdist(non_empty_data, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Ensure the number of rows in pcoa_result and filtered_metadata are the same
if (nrow(pcoa_result$points) == nrow(filtered_metadata)) {
  # Creating a dataframe with PCoA results and corresponding cancer type
  pcoa_coords <- data.frame(PC1 = pcoa_result$points[,1], 
                            PC2 = pcoa_result$points[,2], 
                            CancerType = filtered_metadata$primary_site)
} else {
  stop("The number of rows in PCoA results does not match the number of rows in the filtered metadata.")
}

# Plot PCoA
library(ggplot2)
ggplot(pcoa_coords, aes(x = PC1, y = PC2, color = CancerType)) +
  geom_point(size = 2) +
  labs(title = "PCoA of Cancer Type Genera", x = "PCoA1", y = "PCoA2") +
  theme_minimal()

########################
#Performing tSNE
#Reading the data into R
CommunityMatrix <- read.csv("~/Documents/Input/TabS8.csv")
Metadata <- read.csv("~/Documents/Input/Metadata-TCGA-All-18116-Samples.csv")
#loading dependent package
library(dplyr)

#Removing contaminant genera from CommunityMatrix
#List of genera to remove
ReductGenera <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoeflea", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftia", "Duganella", "Herbaspirillum",                "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobacter", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "Aeromicrobium", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Paenibacillus",                   "Streptococcus", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Pedobacter", "Wautersiella", "Deinococcus")
#Iterating through ReductedGenera to remove contaminant genera(ReductGenera) if they exist
for (genus in ReductGenera) {
  if (genus %in% colnames(CommunityMatrix)) {
    CommunityMatrix <- CommunityMatrix %>% select(-intersect(ReductGenera, colnames(CommunityMatrix)))
  }
}

#Applying read threshold of 10 and any value <10 to 0, any NA to 0
CommunityMatrix <- CommunityMatrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | . < 10, 0, .)))

# Filtering metadata to contain only primary tumor samples
colnames(Metadata)[1] <- "KnightID"
filteredMetadata <- Metadata %>%
  filter(sample_type == 'Primary Tumor') %>%
  select(KnightID, primary_site, data_submitting_center_label)

# Sorting Metadata respective to cancer type in 'primary_site'
filteredMetadata <- filteredMetadata %>%
  mutate(case = case_when(primary_site %in% CancerTypes ~ primary_site, TRUE ~ 'Other'))

# Reducing metadata to retain only KnightID, primary_site and "data_submitting_center_label" 
columns_to_retain <- c("KnightID", "primary_site", "data_submitting_center_label" )
ReducedMetadata <- filteredMetadata %>%
  select(all_of(columns_to_retain))

#Merging CommunityMatrix with Metadata by KnightID the metadata with the raw community matrix
combined_data2 <- merge(CommunityMatrix, ReducedMetadata, by = "KnightID")
colnames(combined_data2)[996] <- "data_submitting_center"

#loading dependent packages
if (!requireNamespace("ClassDiscovery", quietly = TRUE)) install.packages("ClassDiscovery")
if (!requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")
library(ClassDiscovery)
library(Rtsne)

# Ensuring 'genus_data' contains only numerical data
genus_data <- as.matrix(genus_data)

# Removing columns with zero variance
genus_data <- genus_data[, apply(genus_data, 2, var) != 0]

# Verifying that genus_data is not empty after removal
if(ncol(genus_data) == 0) {
  stop("No variance in the genus data; all columns removed.")
}

# Computing Spearman correlation distance matrix
correlation_matrix <- cor(t(genus_data), method = "spearman", use = "pairwise.complete.obs")

# Replacing any remaining NAs in the correlation matrix
correlation_matrix[is.na(correlation_matrix)] <- 0

# Converting to distance matrix
spearman_dist <- as.dist(1 - correlation_matrix)

# Checking for any remaining NAs or non-finite values
if(any(!is.finite(spearman_dist))) {
  stop("The distance matrix contains non-finite values. Please check your data.")
}

# Seting random seed for reproducibility
set.seed(8)

# Performing t-SNE
tsne_result <- Rtsne(spearman_dist, is_distance = TRUE, perplexity = 30, verbose = TRUE)

# Extracting t-SNE coordinates
tsne_coords <- tsne_result$Y

# Ensuring 'combined_data2' has the same number of rows as 'genus_data'
if(nrow(genus_data) != nrow(combined_data)) {
  stop("Mismatch in number of rows between 'genus_data' and 'combined_data'.")
}

# Prepare a data frame for plotting
tsne_df <- data.frame(
  tSNE1 = tsne_coords[,1],
  tSNE2 = tsne_coords[,2],
  TumorType = combined_data2$primary_site,
  SubmittingCenter = combined_data2$data_submitting_center
)

# Plotting t-SNE colored by Tumor Type
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = TumorType)) +
  geom_point(size = 2) +
  labs(title = "t-SNE Plot by Tumor Type", x = "t-SNE1", y = "t-SNE2") +
  theme_minimal()

# Plotting t-SNE colored by Submitting Center
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = SubmittingCenter)) +
  geom_point(size = 2) +
  labs(title = "t-SNE Plot by Submitting Center", x = "t-SNE1", y = "t-SNE2") +
  theme_minimal()


# Filter the data for Thyroid and Lung samples for further PCA as they have maximum PCA variance
thyroid_lung_data <- combined_data %>%
  filter(primary_site %in% c("Thyroid", "Lung"))

# Selecting all numeric columns
numeric_data <- thyroid_lung_data %>% select_if(is.numeric)

# Remove any constant columns (columns with zero variance)
non_constant_data <- numeric_data[, apply(numeric_data, 2, var, na.rm = TRUE) > 0]

# Standardizing data
standardised_data <- scale(non_constant_data)

# Performing PCA
pca_result <- prcomp(standardised_data, center = TRUE, scale. = TRUE)

# Create a dataframe for PCA plotting
pca_df <- data.frame(pca_result$x, 
                     KnightID = thyroid_lung_data$KnightID, 
                     PrimarySite = thyroid_lung_data$primary_site)

# Find samples with the max and min values for PC1 and PC2
max_PC1_sample <- pca_df[which.max(pca_df$PC1), ]
min_PC1_sample <- pca_df[which.min(pca_df$PC1), ]
max_PC2_sample <- pca_df[which.max(pca_df$PC2), ]
min_PC2_sample <- pca_df[which.min(pca_df$PC2), ]

# Print the samples with the max and min PC1 and PC2 values
cat("Max PC1 sample:\n")
print(max_PC1_sample)

cat("Min PC1 sample:\n")
print(min_PC1_sample)

cat("Max PC2 sample:\n")
print(max_PC2_sample)

cat("Min PC2 sample:\n")
print(min_PC2_sample)

# Plotting the PCA results with labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = PrimarySite, label = KnightID)) +
  geom_point(size = 2) +
  labs(title = "PCA of Thyroid and Lung Samples", 
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) +
  theme_minimal() +
  geom_text(aes(label = KnightID), hjust = 0, vjust = 1.5, size = 3)
