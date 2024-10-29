# list.files("data") - if returns character(0), check and set the correct path to working directory
# Read and convert the unfiltered_gene_counts.tsv into a dataframe
unfiltered_data <- data.frame(read.delim("data/unfiltered_gene_counts.tsv",sep="\t"))

# Setting Gene_ID as the row name and delete the Gene_ID column from the df
rownames(unfiltered_data) <- unfiltered_data[,1]
unfiltered_data[,1] <- NULL
any(is.na(unfiltered_data)) # returns false - no NA values in unfiltered_data
# Save the original column names
original_colnames <- colnames(unfiltered_data)

# Normalize data by geometric mean, handling zero counts appropriately
normalized_data <- t(apply(unfiltered_data, 1, function(x) {
  # Check if all values are zero
  if (all(x == 0)) {
    return(rep(0, length(x)))  # Return 0 for this row if all are zero
  }
  
  # Calculate the geometric mean for positive values only
  positive_values <- x[x > 0]
  
  # If there are no positive values, return 0 (or any value you deem appropriate)
  if (length(positive_values) == 0) {
    return(rep(0, length(x)))  # Return 0 if no positive values exist
  }
  
  # Calculate the geometric mean of positive values
  gm <- exp(mean(log(positive_values), na.rm = TRUE))
  
  # Normalize the row by the geometric mean
  return(x / gm)  # Normalize the row
}))

# Assign the original column names back to normalized_data
colnames(normalized_data) <- original_colnames

# Check for NA values in normalized_data
any(is.na(normalized_data))  # return FALSE now - no NA present

# Calculate the size factor for each sample as the median of ratios
normalization_factors <- apply(normalized_data, 2, median)
print(normalization_factors) # all zeros, maybe because there are too many 0 gene counts present?

# The below code lead to some NA values in normalized_data
# # Normalize the gene counts by dividing each element of a row by the geometric mean
# normalized_data <- t(apply(unfiltered_data, 1, function(x) {
#   gm <- exp(mean(log(x[x > 0])))  # Calculate geometric mean for positive values
#   x / gm  # Normalize by geometric mean
# }))
# any(is.na(normalized_data)) # returns true - NA present in normalized_data
