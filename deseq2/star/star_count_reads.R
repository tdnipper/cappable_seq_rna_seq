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
View(samples)

# read in star genecounts as matrix (not dataframe)
dir <- normalizePath("/home/ubuntu/blockvolume/cappable_seq_rna_seq/deseq2/star")
infile_path <- file.path(dir, "deseq2_input.csv")
infile <- as.matrix(read.csv(infile_path, row.names="gene_id"))

# drop x column (index from old dataframe)
infile <- infile[, -which(colnames(infile) %in% "X")]
View(infile)


dds_infection <- DESeqDataSetFromMatrix(countData = infile, colData = samples, ~ infection)
dds_condition <- DESeqDataSetFromMatrix(countData = infile, colData = samples, ~ condition)

# Pre filter 
smallest_group_size <- 2
keep_infection <- rowSums(counts(dds_infection) >= smallest_group_size) >= 10
dds_infection <- dds_infection[keep_infection, ]

#Set up factors
dds_infection$infection <- relevel(dds_infection$infection, ref = "mock")

# Calculate differential expression
dds_infection <- DESeq(dds_infection)
res_infection <- results(dds_infection)
summary(res_infection)

# Filter by padj values to get significant results
res05_infection <- results(dds_infection, alpha = 0.05)
sum(res05_infection$padj < 0.05, na.rm = TRUE) # Returns 2425 genes
summary(res05_infection)

# Export for plotting
res05_infection <- res05_infection[order(-res05_infection$log2FoldChange), ]
write.csv(res05_infection, file = "deseq2/star/geneCounts_infection.csv", row.names = TRUE)
