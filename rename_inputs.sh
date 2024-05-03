#! /bin/bash

# Rename files in raw_data directory to remove the sequencing core processing info from the file name


for file in $(find raw_data -name ".fastq.gz"); do
    basename=$(basename $file)
    newname=${basename%%_*}
    new_path=$(dirname $file)/${newname}.fastq.gz
    mv $file $new_path
done