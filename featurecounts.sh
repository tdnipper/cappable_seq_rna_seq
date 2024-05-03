#! /bin/bash

sudo mkdir -p featurecounts
sudo chmod -R 777 featurecounts

for file in $(find ribodepleted_star_alignments/ -name "*_ribodepleted_sorted.bam"); do
    g=$(basename "$file" _ribodepleted_sorted.bam)
    echo "Counting $p"
    featureCounts -a genome/annotation.gtf \
        -o featurecounts/${g}_ribodepleted_featurecounts.txt \
        -T 4 \
        $file \
        -f \
        -a genome/genomic.gtf \
        -t exon,CDS,gene,transcript
done

for file in $(find star_alignments/ -name "*_sorted.bam"); do
    g=$(basename "$file" _sorted.bam)
    echo "Counting $p"
    featureCounts -a genome/annotation.gtf \
        -o featurecounts/${g}_featurecounts.txt \
        -T 4 \
        $file \
        -f \
        -a genome/genomic.gtf \
        -t exon,CDS,gene,transcript
done