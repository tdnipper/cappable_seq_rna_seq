library("DESeq2")
library("tools")

# Needs samples.txt with sample names in $names column
samples <- read.csv("samples.csv", header = TRUE)
samples <- data.frame(samples)
samples$names <- as.character(samples$names)

# sample$names in alphabetical order to match order of columns from infile
samples <- samples[order(samples$names), ]
samples$condition <- factor(samples$condition)
samples$infection <- factor(samples$infection)
samples$group <- paste(samples$infection, samples$condition, sep = "_")
samples$group <- factor(samples$group)
View(samples)
# read in star genecounts as matrix (not dataframe)
dir <- normalizePath("/home/ubuntu/blockvolume/cappable_seq_rna_seq/deseq2/star")
infile_path <- file.path(dir, "deseq2_input.csv")
infile <- as.matrix(read.csv(infile_path, row.names="gene_id"))

# drop x column (index from old dataframe)
infile <- infile[, -which(colnames(infile) %in% "X")]

# Set up DE experiment using infection status as the main condition
# Use cappable-seq enrichment as additional factor
# dds_infection <- DESeqDataSetFromMatrix(countData = infile, colData = samples, ~ condition + infection)
dds_infection <- DESeqDataSetFromMatrix(countData = infile, colData = samples, ~ group)

# Specify the reference levels in each experimental factor
# dds_infection$infection <- relevel(dds_infection$infection, ref = "mock")
# dds_infection$condition <- relevel(dds_infection$condition, ref = "input")

# Pre filter genes that have < 10 reads across 2 samples
smallest_group_size <- 2
keep_infection <- rowSums(counts(dds_infection) >= smallest_group_size) >= 10
dds_infection <- dds_infection[keep_infection, ]

# Calculate differential expression
dds_infection <- DESeq(dds_infection)
resultsNames(dds_infection)
res_infection <- results(dds_infection, contrast = c("group", "mock_enriched", "infected_enriched"))

# Apply FDR filter of 0.05 to results 
res05_infection <- results(dds_infection, alpha = 0.05)
sum(res05_infection$padj < 0.05, na.rm = TRUE) # Returns # of genes that pass FDR alpha
summary(res05_infection)

# Export for plotting in csv
res05_infection <- res05_infection[order(-res05_infection$log2FoldChange), ]
write.csv(res05_infection, file = "deseq2/star/geneCounts_infection.csv", row.names = TRUE)
significant <- subset(res05_infection, padj <= 0.05)
write.csv(significant, file = "deseq2/star/geneCounts_infection_significant.csv", row.names = TRUE)
