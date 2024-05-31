library("DESeq2")
library("tximeta")
library("tools")

# Needs samples.txt with sample names in $names column
samples <- read.csv("samples.csv", header = TRUE)
samples <- data.frame(samples)
samples$condition <- factor(samples$condition)
samples$infection <- factor(samples$infection)

# read in star genecounts
dir <- normalizePath("/home/ubuntu/blockvolume/cappable_seq_rna_seq/deseq2_star")
infile <- file.path(dir, "deseq2_input.csv")

dds_infection <- DESeqDataSetFromMatrix(countData = infile, colData = samples, ~ infection)
