if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
install.packages("pheatmap")


library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library("pheatmap")
library(ggrepel)
library(ggplot2)

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
integer_data <- round(filtered_counts)
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
res_df <- as.data.frame(res) %>% dplyr::mutate(significant = padj < 0.5)
res_df$Gene <- rownames(res_df)
res_labelled <- res_df %>%
  dplyr::filter(significant)

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
plotMA(res, ylim=c(-10, 25))
# using the shrunken log2 fold changes for our plot
plotMA(resLFC, ylim=c(-24, 6))

ggplot(res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
  geom_point() +
  ggrepel::geom_label_repel(data=res_labelled, aes(label=Gene), max.overlaps=15) +
  #geom_text(aes(label = Gene))
  theme_bw()

# Convert the LFC DESeq2 results to a dataframe
res_lfc_df <- as.data.frame(resLFC) %>% dplyr::mutate(significant = padj < 0.5)
res_lfc_df$Gene <- rownames(res_lfc_df)
res_lfc_labelled <- res_lfc_df %>%
  dplyr::filter(significant)

ggplot(res_lfc_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
  geom_point() +
  ggrepel::geom_label_repel(data=res_lfc_labelled, aes(label=Gene), max.overlaps=15) +
  #geom_text(aes(label = Gene))
  theme_bw()


#heatmap of most significant gene counts
significant <- rownames(na.omit(res)[na.omit(res)$padj <= 0.0005,])
pheatmap(normalized_counts[rownames(normalized_counts) %in% significant,])
# TODO use z scores


# get summary of results
summary(res)

# see results sorted by adjusted p value
res <- res[order(res$padj),]
head(res)

# the counts of reads for a single gene across the groups
ENSG00000223019 <- plotCounts(analysis, gene="ENSG00000223019", intgroup="Sex", 
                returnData=TRUE)

#violin plot for ENSG00000223019
ggplot(ENSG00000223019, aes(x=Sex, y=count), ) +
  geom_violin(trim=FALSE) + 
  stat_summary(fun.y=mean, geom="point", size=4, aes(shape="mean")) + 
  stat_summary(fun.y=median, geom="point", size=10, color="red", aes(shape="median")) + 
  scale_shape_manual("", values=c("mean"="x", "median"="•")) +
  ggtitle("ENSG00000223019 gene counts by sex")


# the counts of reads for a single gene across the groups
ENSG00000252353 <- plotCounts(analysis, gene="ENSG00000252353", intgroup="Sex", 
                              returnData=TRUE)

#violin plot for ENSG00000252353
ggplot(ENSG00000252353, aes(x=Sex, y=count), ) +
  geom_violin(trim=FALSE) + 
  stat_summary(fun.y=mean, geom="point", size=4, aes(shape="mean")) + 
  stat_summary(fun.y=median, geom="point", size=10, color="red", aes(shape="median")) + 
  scale_shape_manual("", values=c("mean"="x", "median"="•")) +
  ggtitle("ENSG00000252353 gene counts by sex")



# the counts of reads for a single gene across the groups
MI0016797 <- plotCounts(analysis, gene="MI0016797", intgroup="Sex", 
                              returnData=TRUE)
MI0016797
#violin plot for MI0016797
ggplot(MI0016797, aes(x=Sex, y=count)) +
  geom_violin(scale="width", trim=FALSE) + 
  stat_summary(fun.y=mean, geom="point", size=4, aes(shape="mean")) + 
  stat_summary(fun.y=median, geom="point", size=10, color="red", aes(shape="median")) + 
  scale_shape_manual("", values=c("mean"="x", "median"="•")) +
  ggtitle("MI0016797 gene counts by sex")


ENSG00000194297 <- plotCounts(analysis, gene="ENSG00000194297", intgroup="Sex", 
                              returnData=TRUE)

#violin plot for ENSG00000194297
ggplot(ENSG00000194297, aes(x=Sex, y=count), ) +
  geom_violin(scale="width", trim=FALSE) + 
  stat_summary(fun.y=mean, geom="point", size=4, aes(shape="mean")) + 
  stat_summary(fun.y=median, geom="point", size=10, color="red", aes(shape="median")) + 
  scale_shape_manual("", values=c("mean"="x", "median"="•")) +
  ggtitle("ENSG00000194297 gene counts by sex")


MIMAT0000757 <- plotCounts(analysis, gene="MIMAT0000757", intgroup="Sex", 
                              returnData=TRUE)

#violin plot for MIMAT0000757
ggplot(MIMAT0000757, aes(x=Sex, y=count), ) +
  geom_violin(scale = "width", trim=FALSE) + 
  stat_summary(fun.y=mean, geom="point", size=4, aes(shape="mean")) + 
  stat_summary(fun.y=median, geom="point", size=10, color="red", aes(shape="median")) + 
  scale_shape_manual("", values=c("mean"="x", "median"="•")) +
  ggtitle("MIMAT0000757 gene counts by sex")


# ENSG00000223019
# ENSG00000252353
# MI0016797
# ENSG00000194297
# MIMAT0000757


mcols(res)$description


