unfiltered_data <- data.frame(read.delim("data/unfiltered_gene_counts.tsv",sep="\t"))
rownames(unfiltered_data) <- unfiltered_data[,1]
unfiltered_data[,1] <- NULL
head(unfiltered_data)

normalized_data <- t(apply(unfiltered_data, 1, function(x) x / mean(x)))
head(normalized_data)
