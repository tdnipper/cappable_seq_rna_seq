#! /bin/bash

# Run salmon on STAR alignments in alignment-based mode
# Star alignments must be aligned to transcriptome
# Salmon doesn't need an index here

source env.sh

if [ ! -d ribodepleted_star_alignments ]; then
    echo "Run STAR first!"
    exit 1
fi

if [ ! -d ribodepelted_star_alignments/salmon_quantifications ]; then
    mkdir ribodepleted_star_alignments/salmon_quantifications
    sudo chown $USER_ID:$GROUP_ID ribodepleted_star_alignments/salmon_quantifications
fi

for file in $(find ribodepleted_star_alignments -name "Aligned.toTranscriptome.out.bam"); do
    name=$(basename $file "ribodepleted_Aligned.toTranscriptome.out.bam")
    echo "Quantifying ${name}"
    salmon quant \
\        -l A \
        -r $file \
        --validateMappings \
        -p 4 \
        -o ribodepleted_star_alignments/salmon_quantifications/${name}
done

# Rename files with sample name
for file in $(find ribodepleted_star_alignments/salmon_quantifications -name "quant.sf"); do
    name=$(dirname $file | awk -F/ '{print $NF}')
    mv $file ribodepleted_star_alignments/salmon_quantifications/${name}_quant.sf
done