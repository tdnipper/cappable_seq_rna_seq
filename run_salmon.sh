#! /bin/bash

source env.sh

if [ ! -d salmon_quantifications ]; then
    mkdir salmon_quantifications
    sudo chown $USER_ID:$USER_ID salmon_quantifications
fi

if [ ! -d genome/salmon_index ]; then
salmon index \
    -t genome/hybrid_gentrome.fasta \
    -d genome/decoys.txt \
    -p 4 \
    -i genome/salmon_index
fi

for file in $(find ribodepleted_reads/ -name "*trimmed_nonrRNA.fq.gz"); do
    name=$(basename $file "_trimmed_nonrRNA.fq.gz")
    echo "Quantifying ${name}"
    salmon quant \
        -i genome/salmon_index \
        -l A \
        -r $file \
        --validateMappings \
        -p 4 \
        -o salmon_quantifications/${name}
done

