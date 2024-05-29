# install.packages("BiocManager")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("DESeq2", dependencies = TRUE)
# BiocManager::install("tximport", dependencies = TRUE)
# BiocManager::install("tximeta", dependencies = TRUE)
# BiocManager::install("BiocFileCache", dependencies = TRUE)
# Start here once packages are installed

library(DESeq2)
library(tximeta)
library(tools)
library(BiocFileCache)

# Needs samples.txt with sample names in $names column
samples <- read.csv("samples.csv", header = TRUE)
samples <- data.frame(samples)

# Needs salmon_quantification directory with salmon quant.sf files
dir <- normalizePath("salmon_quantification")
files <- file.path(dir, paste0(samples$names, "/", samples$names, "_quant.sf"))
samples$files <- files
file.exists(samples$files)

# Make custom linked txome from salmon index to get hybrid WSN/human

salmon_dir <- normalizePath("genome/salmon_index")
transcript_dir <- normalizePath("genome/hybrid_exon.fasta")
genome_dir <- normalizePath("genome/hybrid_genome.fasta")
gtf_dir <- normalizePath("genome/hybrid_annotated_gffread.gtf")

txome <- makeLinkedTxome(
  indexDir=salmon_dir,
  source="GENCODE homemade",
  organism="human/WSN",
  release="GRCh38.v46",
  genome=genome_dir,
  fasta=transcript_dir,
  gtf=gtf_dir,
  write=TRUE,
  jsonFile="linkedTxomeTbl.json"
)
# bfcloc <- getTximetaBFC()
# bfc <- BiocFileCache(bfcloc)
# bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)

# Make a summarized experiment object using provided metadata and salmon quant
se <- tximeta(samples)
colData(se)
assayNames(se)
rowRanges(se)
seqinfo(se)

# Get transcript db
edb <- retrieveDb(se)
class(edb)

# Get exons
# se_exons <- addExons(se)
# rowRanges(se_exons[[1]])
#### This doesn't work for some reason, but doesn't seem to affect gse output

# Summarize to gene level for subsequent quantification
gse <- summarizeToGene(se)
rowRanges(gse)
dim(gse)

# Import to Deseq2
dds_condition <- DESeqDataSet(gse, design = ~ condition)
dds_infection <- DESeqDataSet(gse, design = ~ infection)

# Pre-filter rows with < 10 reads in 2 samples
smallest_group_size <- 2
keep_condition <- rowSums(counts(dds_condition) >= smallest_group_size) >= 10
dds_condition <- dds_condition[keep_condition,]
keep_infection <- rowSums(counts(dds_infection) >= smallest_group_size) >= 10
dds_infection <- dds_infection[keep_infection,]

# Set up factors
levels(dds_condition$condition)
dds_condition$condition <- relevel(dds_condition$condition, ref = "input")
dds_infection$infection <- relevel(dds_infection$infection, ref = "mock")
# Calculate differential expression
dds_condition <- DESeq(dds_condition)
res_condition <- results(dds_condition)
# resultsNames(dds_condition)
resLFC_condition <- lfcShrink(dds_condition, coef = "condition_enriched_vs_input", type = "apeglm")
# res_condition
# summary(res_condition)
dds_infection <- DESeq(dds_infection)
# resultsNames(dds_infection)
resLFC_infection <- lfcShrink(dds_infection, coef = "infection_infected_vs_mock", type = "apeglm")
res_infection <- results(dds_infection)
# summary(res_infection)

# Filter padj values to get significant results
res05_condition <- results(dds_condition, alpha = 0.05)
sum(res05_condition$padj < 0.05, na.rm = TRUE)
summary(res05_condition)
res05_infection <- results(dds_infection, alpha = 0.05)
sum(res05_infection$padj < 0.05, na.rm = TRUE)
summary(res05_infection)

# Plot
res05_condition <- res05_condition[order(-res05_condition$log2FoldChange),]
res05_infection <- res05_infection[order(-res05_infection$log2FoldChange),]

# Export res05_infection to CSV
write.csv(res05_condition, file = "res05_condition.csv", row.names = TRUE)
write.csv(res05_infection, file = "res05_infection.csv", row.names = TRUE)
