install.packages("BiocManager")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("tximport", dependencies = TRUE)
BiocManager::install("tximeta", dependencies = TRUE)

# Start here once packages are installed

library(DESeq2)
library(tximeta)
library(tools)

# Needs samples.txt with sample names in $names column
samples <- read.table("samples.txt", header = TRUE)
samples <- data.frame(samples)

# Needs salmon_quantification directory with salmon quant.sf files
dir <- normalizePath("salmon_quantification")
files <- file.path(dir, paste0(samples$names, "/", samples$names, "_quant.sf"))
samples$files <- files
View(samples)
file.exists(samples$files)

# Make custom linked txome from salmon index to get hybrid WSN/human
txome <- makeLinkedTxome(
    "genome/salmon_index",
    "de-novo",
    "Homo sapiens and WSN",
    "GENCODE/custom", "GRCh38/WSN",
    "genome/hybrid_exon.fasta",
    "genome/hybrid_annotated_cat.gtf",
    "tx.json"
)
makeLinkedTxome(txome)

se <- tximeta(samples)
