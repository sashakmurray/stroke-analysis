if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")

library(dplyr)
library(DESeq2)
library(EnhancedVolcano)


# list.files("data") - if returns character(0), check and set the correct path to working directory
# Read and convert the unfiltered_gene_counts.tsv into a dataframe
unfiltered_data <- data.frame(read.delim("data/unfiltered_gene_counts.tsv",sep="\t", row.names=1))

# filter out lincRNAs
rna_biotypes <- read.csv('data/rna_biotypes.csv', header = TRUE, row.names=1)
non_lincRNA <- rownames(rna_biotypes[rna_biotypes$Biotype != "lincRNA", , drop = FALSE])
counts_without_lincRNA <- unfiltered_data[rownames(unfiltered_data) %in% non_lincRNA,]

#Remove genes where counts are <10 across all samples
filtered_counts <- counts_without_lincRNA[rowSums(counts_without_lincRNA >= 10) > 0, ]

average_counts <- rowMeans(filtered_counts)
cutoff <- quantile(average_counts, 0.15)
filtered_counts <- filtered_counts[average_counts > cutoff, ]

# convert gene counts to integers
integer_data <- round(unfiltered_data)
head(integer_data)

# use annotations file to get metadata for analysis of gene counts
metaData <- read.csv('data/annotations.tsv', header = TRUE, sep = "\t", row.names=1)
head(metaData)

metaData$Sex <- as.factor(metaData$Sex)



# using DESeq2 to generate size factors
dds <- DESeqDataSetFromMatrix(countData = integer_data, colData = metaData, tidy = FALSE, design = ~Sex)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
# write normalized data to csv
write.table(normalized_counts, file="data/normalized_counts.csv", sep=",")

# run DESeq analysis
analysis <- DESeq(dds)

# Extract results for Male vs. Female comparison
res <- results(analysis, contrast = c("Sex", "0", "1"))
# Add significance for visualization
res$significance <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")

# these are Log fold change shrinkage results
# removes the noise associated with log2 fold changes from 
# low count genes without requiring arbitrary filtering thresholds.
resLFC <- lfcShrink(analysis, coef="Sex_1_vs_0", type="ashr")

# Convert the DESeq2 results to a dataframe
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

# Create the volcano plot using EnhancedVolcano
EnhancedVolcano(res,
                lab = rownames(res),  # Labels for the genes
                x = 'log2FoldChange',  # Column for x-axis
                y = 'padj',            # Column for y-axis (adjusted p-value)
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 Adjusted p-value',
                title = 'Volcano Plot of DEGs between Male and Female',
                pCutoff = 0.05,       # Significance cutoff
                FCcutoff = 0,         # Fold change cutoff
                pointSize = 3.0,      # Size of points
                labSize = 3.0,        # Size of labels
                gridlines.major = FALSE, # Disable major gridlines
                gridlines.minor = FALSE  # Disable minor gridlines
)


# MA-plot shows the log2 fold changes attributable to a given 
#variable over the mean of normalized counts for all the samples
plotMA(res, ylim=c(-7, 7))
# using the shrunken log2 fold changes for our plot
plotMA(resLFC, ylim=c(-6, 6))


