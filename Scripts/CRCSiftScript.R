#Reading the data into R
CommunityMatrix <- read.csv("~/Documents/Assignments/r_work/Assignemnt1Done/TabS8.csv")
Metadata <- read.csv("~/Documents/Assignments/r_work/Assignemnt1Done/Metadata-TCGA-All-18116-Samples.csv")
#loading dependent package
library(dplyr)

#Removing contaminant genera from CommunityMatrix
#List of genera to remove
ReductGenera <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoeflea", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftia", "Duganella", "Herbaspirillum",                "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobacter", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "Aeromicrobium", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Paenibacillus",                   "Streptococcus", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Pedobacter", "Wautersiella", "Deinococcus")

#Iterating through ReductedGenera to remove contaminant genera if they exist
for (genus in ReductGenera) {
  if (genus %in% colnames(CommunityMatrix)) {
    CommunityMatrix <- CommunityMatrix %>% select(-all_of(genus))
  }
}

#Applying read threshold of 10 and any value <10 to 0, any NA to 0
CommunityMatrix <- CommunityMatrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | . < 10, 0, .)))

#Filtering metadata to contain only primary tumour samples
colnames(Metadata)[1] <- "KnightID"
filteredMetadata <- Metadata %>% filter(sample_type=='Primary Tumor') %>% select(KnightID, primary_site) 
filteredMetadata <- filteredMetadata %>% mutate(crc = case_when( primary_site == 'Colorectal' ~ 'Colorectal', TRUE ~ 'Other' )) 
colnames(filteredMetadata)[2] <- "primary_site"
filteredMetadata <- filteredMetadata %>% select(-primary_site) 

#Deleting all other FilteredMetadata columns except KnightID and primary_site
colnames(filteredMetadata)[1] <- "KnightID"
columns_to_retain <- c("KnightID", "primary_site")
colnames(filteredMetadata)[2] <- "primary_site"
ReducedMetadata <- filteredMetadata %>%
  select(all_of(columns_to_retain))

#Merging CommunityMatrix with Metadata by KnightID the metadata with the raw community matrix
combined_data <- merge(CommunityMatrix, ReducedMetadata, by = "KnightID")


#seperating combined_data into Colorectal and all other cancer types
Colorectal <- "Colorectal"
# Dataset containing information about Colorectal
Colorectaldata <- combined_data %>%
  filter(primary_site == Colorectal)
# Dataset containing information about all other cancer types
AllOtherCancerTypes <- combined_data %>%
  filter(primary_site != Colorectal)

#looping through each genera to calc no of samples >0 and those  equal to zero for colorectal 
results <- data.frame(Genus = character(), crc_pos = integer(), crc_neg = integer(), stringsAsFactors = FALSE)

for (genus in colnames(Colorectaldata)) {
  crc_pos <- sum(Colorectaldata[[genus]] > 0, na.rm = TRUE)  # Number of samples > 0
  crc_neg <- sum(Colorectaldata[[genus]] == 0, na.rm = TRUE) # Number of samples == 0
  
  # Adding the results to the data frame
  results <- rbind(results, data.frame(Genus = genus, crc_pos = crc_pos, crc_neg = crc_neg))
}

#looping through each genera to calc no of samples >0 and those  equal to zero for othercancertypes
result2 <- data.frame(Genus = character(), other_pos = integer(), other_neg = integer(), stringsAsFactors = FALSE)

for (genus in colnames(AllOtherCancerTypes)) {
  other_pos <- sum(AllOtherCancerTypes[[genus]] > 0, na.rm = TRUE)  # Number of samples > 0
  other_neg <- sum(AllOtherCancerTypes[[genus]] == 0, na.rm = TRUE) # Number of samples == 0
  
  # Adding the results to the data frame
  result2 <- rbind(result2, data.frame(Genus = genus, other_pos = other_pos, other_neg = other_neg))
}
# merging results(colorectal) with result2 (allothercancers) by Genus
FinalResults <- merge(results, result2, by = "Genus")
print(FinalResults)
write.csv(FinalResults, "FinalResults.csv", row.names = FALSE)

#Performing Fisher's exact test and return p-value
get_p_value <- function(crc_pos, crc_neg, other_pos, other_neg) {
  # Creating matrix for Fisher's test
  matrix_data <- matrix(c(crc_pos, crc_neg, other_pos, other_neg), nrow = 2, byrow = TRUE)
  
  # Checking matrix dimensions (should be 2x2)
  if (all(dim(matrix_data) == c(2, 2))) {
    test_result <- fisher.test(matrix_data)
    return(test_result$p.value)
  } else {
    # Return NA if the matrix is not 2x2
    return(NA)
  }
}

# Applying function to each row and creating new column for p-values
FinalResults <- FinalResults %>%
  rowwise() %>%
  mutate(p_value = get_p_value(crc_pos, crc_neg, other_pos, other_neg)) %>%
  ungroup()
#adjusting p_values
p.adjust(FinalResults$p_value, method='BH')
# Sorting dataframe by p-value in descending order
sorted_FinalResults <- FinalResults %>%
  arrange(p_value)

# ViewS sorted dataframe
print(sorted_FinalResults)
write.csv(sorted_FinalResults, "DecontaminatedSorted_FinalResults.csv", row.names = FALSE)



