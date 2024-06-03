library("DESeq2")
library("tximeta")
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
