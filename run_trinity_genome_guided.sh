#! /bin/bash

# genome guided trinity on merged bam file

Trinity --genome_guided_bam trinity/merged_sorted.bam \
    --genome_guided_max_intron 1000000 \
    --max_memory 30G \
    --CPU 4 \
    --output trinity/