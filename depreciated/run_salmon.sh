#! /bin/bash

#####################
# DEPRECIATED IN FAVOR OF find_salmon_replicates.py
#####################

# Run salmon in mapping mode on trimmed, ribodepleted reads
# Circumvents star and gives counts per transcript to pass to tximport directly

source env.sh

if [ ! -d salmon_quantifications ]; then
    mkdir salmon_quantifications
    sudo chown $USER_ID:$GROUP_ID salmon_quantifications
fi

if [ ! -d genome/salmon_index ]; then
salmon index \
    -t genome/hybrid_gentrome.fasta \
    -d genome/decoys.txt \
    -p 4 \
    -i genome/salmon_index
fi

for file in $(find ribodepleted_reads/ -name "*_nonrRNA.fq.gz"); do
    name=$(basename $file "_nonrRNA.fq.gz")
    echo "Quantifying ${name}"
    salmon quant \
        -i genome/salmon_index \
        -l A \
        -r $file \
        --validateMappings \
        -p 4 \
        -o salmon_quantifications/${name} \
        --writeUnmappedNames \
        --seqBias \
        --gcBias \
        --reduceGCMemory
done

# Rename quant.sf files to include sample name from subdir
for file in $(find salmon_quantifications -name "quant.sf"); do
    name=$(dirname $file | awk -F/ '{print $NF}')
    mv $file $(dirname $file)/${name}_quant.sf
done