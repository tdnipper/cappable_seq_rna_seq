# Validate trinity transcriptome against hybrid using bowtie2 alignment

if [! -d "genome/bowtie2_index"] then
    echo("Please make a bowtie index at genome/bowtie2_index first")
    exit 1
fi

bowtie2 -p 4 \
    -x genome/bowtie2_index \
    -U trinity/Trinity-GG.fasta \
    -S trinity/Trinity-GG.bam
    