#Get significant SNPs
# List files matching the new pattern Trait1.assoc.txt, Trait2.assoc.txt, ... Trait27.assoc.txt ; Generated from GEMMA
file_list <- list.files(pattern = "^Trait[0-9]+\\.assoc\\.txt$", full.names = TRUE)
print(file_list)  # check that the correct files are loaded

for (file_path in file_list) {
  original_data <- fread(file_path)
  
  # check
  print(paste("Processing file:", file_path))
  print(head(original_data))
  
  # p value is smaller than X (say 5.02E-8)
  filtered_data <- original_data[original_data[[13]] < 5.02E-8, ]
  
  # extract sig snps --- rows
  if (nrow(filtered_data) > 0) {
    # Extract trait code
    trait_name <- gsub("^sort_|\\.txt$", "", basename(file_path))
    assign(paste0("Sig_", trait_name), filtered_data)
    
    # check trait name correctly assigned
    print(paste("Trait name:", trait_name))
    print(head(filtered_data))
  } else {
    print(paste("No significant data so useless file:", file_path))
  }
}

# Get a list of sig signals names starting with "Sig_"
data_frames <- ls(pattern = "^Sig_")
for (df_name in data_frames) {
  df <- get(df_name)
  
  # Add the "Trait" column to the data frame so i know the corresponding traits
  df$Trait <- gsub("^Sig_", "", df_name)  # Extract the trait name from the data frame name
  
  assign(df_name, df)
  
  # final check
  print(paste("Data frame:", df_name))
  print(head(df))
}

# nope thats the final check :)
print("Data frames created:")
print(data_frames)
df_names <- c("Sig_Trait1.assoc", "Sig_Trait2.assoc", "Sig_Trait3.assoc", "Sig_Trait4.assoc", "Sig_Trait5.assoc", "Sig_Trait6.assoc", "Sig_Trait7.assoc", "Sig_Trait8.assoc", "Sig_Trait9.assoc")
df_list <- list()
# Iterate over the data frame names
for (df_name in df_names) {
  
  # Check if the data frame exists --- significant signals exist :)
  
  if (exists(df_name)) {
    
    # If it exists, retrieve it and add it to the list
    
    df_list[[df_name]] <- get(df_name)
    
  } else {
    
    # If it doesn't exist, FUCK IT
    
    print(paste("Data frame", df_name, "does not exist. Skipping..."))
    
  }
  
}

Sig <- do.call(rbind, df_list)
write.table(Sig, "Sig.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

B <- fread("/data/gpfs/projects/punim1869/users/yunhongh1/workspace/Wild/Filtered_geno0.2_mind0.5_maf0.05_154.blocks.det") #LD blocks across genome generated from plink
formatted_data <- separate_rows(B, SNPS, sep = "\\|") %>%
  select(CHR, BP1, SNPS) #plink block names
colnames(formatted_data) <- c("CHR", "BlockID", "SNP")
write.table(formatted_data, file = "ld_blocks_formatted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

matching_rows <- formatted_data[formatted_data$SNP %in% Sig$rs, ]
matching_rows <- matching_rows %>%
  left_join(Sig[, c("rs", "Trait")], by = c("SNP" = "rs"))  

matching_rows <- matching_rows %>%
  group_by(BlockID) %>%
  summarise(
    CHR = first(CHR),                
    BlockID = first(BlockID),        
    SNP = first(SNP),                
    Traits = paste(unique(Trait), collapse = ";") 
  )

POS <- matching_rows$BlockID
B2 <- B[BP1 %in% POS]
B2_with_traits <- B2 %>%
  left_join(matching_rows[, c("BlockID", "Traits")], by = c("BP1" = "BlockID"))

write.table(B2_with_traits, "Significant_Block_with_Traits.txt", row.names = FALSE, sep = "\t", quote = FALSE)
