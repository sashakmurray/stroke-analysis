BiocManager::install("vsn")
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(vsn)



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
colData(dds)

normalized_counts <- counts(dds, normalized=TRUE)
# write normalized data to csv
write.table(normalized_counts, file="data/normalized_counts.csv", sep=",")

# run DESeq analysis
analysis <- DESeq(dds)

# Extract results for Male vs. Female comparison
res_sex <- results(analysis, contrast = c("Sex", "0", "1"))

# these are Log fold change shrinkage results
# removes the noise associated with log2 fold changes from 
# low count genes without requiring arbitrary filtering thresholds.
resLFC_sex <- lfcShrink(analysis, coef="Sex_1_vs_0", type="ashr")

# Convert the DESeq2 results to a dataframe
res_df <- as.data.frame(res_sex) %>% dplyr::mutate(significant = padj < 0.5)
res_df$Gene <- rownames(res_df)
rna_biotypes$Gene <- rownames(rna_biotypes)
res_df <- merge(res_df, rna_biotypes, by.x = "Gene", by.y = "Gene", all.x = TRUE)
res_labelled <- res_df %>%
  dplyr::filter(significant)


# Creating visualizations below
# Create the volcano plot using EnhancedVolcano: width = 790, height = 625
EnhancedVolcano(res_sex,
                lab = res_df$Biotype,  # Labels for the genes
                x = 'log2FoldChange',  # Column for x-axis
                y = 'padj',            # Column for y-axis (adjusted p-value)
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 Adjusted p-value',
                title = 'Volcano Plot of Differential Expression Genes',
                pCutoff = 0.05,       # Significance cutoff
                FCcutoff = 0,         # Fold change cutoff
                pointSize = 2.5,      # Size of points
                labSize = 3.0,        # Size of labels
                gridlines.major = FALSE, # Disable major gridlines
                gridlines.minor = FALSE  # Disable minor gridlines
)


# MA-plot shows the log2 fold changes attributable to a given 
#variable over the mean of normalized counts for all the samples
plotMA(res, ylim=c(-10, 25))
# using the shrunken log2 fold changes for our plot
plotMA(resLFC, ylim=c(-24, 6))


# labelled MA plot
ggplot(res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
  geom_point() +
  ggrepel::geom_label_repel(data=res_labelled, aes(label=Gene), max.overlaps=15) +
  #geom_text(aes(label = Gene))
  theme_bw()

# Convert the LFC DESeq2 results to a dataframe
res_lfc_df <- as.data.frame(resLFC_sex) %>% dplyr::mutate(significant = padj < 0.5)
res_lfc_df$Gene <- rownames(res_lfc_df)
res_lfc_labelled <- res_lfc_df %>%
  dplyr::filter(significant)

# labelled MA plot with LFC genes
ggplot(res_lfc_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
  geom_point() +
  ggrepel::geom_label_repel(data=res_lfc_labelled, aes(label=Gene), max.overlaps=15) +
  #geom_text(aes(label = Gene))
  theme_bw()


# get summary of results
summary(res_sex)

# see results sorted by adjusted p value
res <- res[order(res_sex$padj),]
head(res_sex)


# Define color mapping for annotations
annotation_colors <- list(Sex = c("0" = "#11D6E6", "1" = "#F25E52"))
# Create annotation for samples
annotate_sex <- as.data.frame(colData(dds)[, "Sex", drop=FALSE])

# Exclude rows where 'padj' is NA, then filter for significant genes
significant_genes <- rownames(res_sex[!is.na(res_sex$padj) & res_sex$padj < 0.05 & abs(res_sex$log2FoldChange) > 1, ])

# Transform the filtered count data for visualization 
vsd <- vst(dds, blind=FALSE)
transformed_vst <- assay(vsd)[significant_genes, ]
head(transformed_vst, 3)
meanSdPlot(assay(vsd))
# Compute row-wise z-scores
z_scores_vst <- t(scale(t(transformed_vst)))
# Remove rows with NA values (e.g., genes with zero variance)
z_scores_vst <- z_scores_vst[complete.cases(z_scores_vst), ]
# Creating heatmaps of the transformed data
z_scores_vst <- as.data.frame(z_scores_vst)
z_scores_vst$Gene <- rownames(z_scores_vst)
z_scores_vst <- merge(z_scores_vst, rna_biotypes, by.x = "Gene", by.y = "Gene", all.x = TRUE)
rownames(annotate_sex) == colnames(z_scores_vst)[!colnames(z_scores_vst) %in% c("Gene", "Biotype")]  # Ensure alignment

# width = 885, height = 672
pheatmap(z_scores_vst[!colnames(z_scores_vst) %in% c("Gene", "Biotype")], 
         cluster_rows=TRUE, 
         labels_row = z_scores_vst$Biotype,
         cluster_cols=FALSE,
         annotation_col=annotate_sex,
         annotation_colors = annotation_colors,
         main = "Heatmap of Significant Genes (FDR < 0.05)")

# Another type of transformation: log2(n + 1)
ntd <- normTransform(dds)
transformed_ntd <- assay(ntd)[significant_genes, ]
meanSdPlot(transformed_ntd)
z_scores_ntd <- t(scale(t(transformed_ntd)))
z_scores_ntd <- z_scores_ntd[complete.cases(z_scores_ntd), ]
z_scores_ntd <- as.data.frame(z_scores_ntd)
z_scores_ntd$Gene <- rownames(z_scores_ntd)
z_scores_ntd <- merge(z_scores_ntd, rna_biotypes, by.x = "Gene", by.y = "Gene", all.x = TRUE)
rownames(annotate_sex) == colnames(z_scores_ntd)[!colnames(z_scores_ntd) %in% c("Gene", "Biotype")]  # Ensure alignment

pheatmap(z_scores_ntd[!colnames(z_scores_ntd) %in% c("Gene", "Biotype")],
         cluster_rows=TRUE, 
         labels_row = z_scores_ntd$Biotype,
         cluster_cols=FALSE,
         annotation_col=annotate_sex,
         annotation_colors = annotation_colors,
         main = "Heatmap of Significant Genes (FDR < 0.05)")

# TODO: weighted gene correlation network
# TODO: box-plot / violin plot for individual genes
# TODO: visualizations for demographics



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

violin plot for MI0016797
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

