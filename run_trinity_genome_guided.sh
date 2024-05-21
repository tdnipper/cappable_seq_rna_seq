#! /bin/bash

CPU=$(nproc)

# genome guided trinity on merged bam file

Trinity --genome_guided_bam trinity_genome_guided/merged_sorted.bam \
    --genome_guided_max_intron 1000000 \
    --max_memory 30G \
    --CPU $CPU \
    --output trinity_genome_guided/