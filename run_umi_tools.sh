#! /bin/bash

source env.sh

# Deduplicate all sorted BAM files in the directory
for file in $(find ribodepleted_star_alignments/ -name "*_sorted.bam"); do
    name=$(basename $file "_sorted.bam")
    echo "Deduplicating $name" 
    umi_tools dedup -I $file -S ${file%.*}_dedup.bam \
        --method=directional \
        --extract-umi-method=read_id \
        --umi-separator=rbc: \
        --log ${file%.*}_dedup.log \
        --error ${file%.*}_dedup_error.log
done

# Index sorted dedup files
for file in $(find ribodepleted_star_alignments/ -name "*_sorted_dedup.bam"); do
    echo "Indexing $file"
    samtools index $file
done
