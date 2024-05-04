#! /bin/bash

# Rename files in raw_data directory to remove the sequencing core processing info from the file name

for file in $(find raw_data -name "*.fastq.gz"); do 
    basename=$(basename $file ".fastq.gz")
    newname=${basename%%_*} # remove everything after first '_'
    replicate=${basename%_*} # remove everything after basename last '_'
    replicate=${replicate##*_} # remove everything before replicate last '_'
    new_path=$(dirname $file)/${newname}_${replicate}.fastq.gz # construct new path
    mv $file $new_path # rename file
done