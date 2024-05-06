#! /bin/bash

source env.sh

if [ ! -d salmon_quantifications ]; then
    mkdir salmon_quantifications
    sudo chown $USER_ID:$USER_ID salmon_quantifications
fi

if [! -f salmon_index ]; then
    salmon index -t genome/transcriptome.fasta \
        -i salmon_quantifications/salmon_index
fi

for file in $(find ribodepleted_star_alignments/ -name "*_transcriptome_sorted.bam"); do
    name = $(basename $file "_ribodepleted_transcriptome_sorted.bam")
    echo "Quantifying ${name}"
    salmon quant \
        -i salmon_quantifications/salmon_index \
        -l A \
        -1 $file \

