#! /bin/bash

source env.sh

sudo mkdir -p $BASEDIR/qc_reports/samtools-stats
sudo chown -R $USER_ID:$GROUP_ID $BASEDIR/qc_reports/samtools-stats
sudo chmod -R 775 /home/ubuntu/blockvolume/dusp11_clip-seq/qc_reports/samtools-stats

for file in $(find $BASEDIR -name "*_sorted.bam"); do
    name=$(basename $file "_sorted.bam")
    echo "Processing $name" 
    samtools stats $file > $BASEDIR/qc_reports/samtools-stats/${name}.stats --threads 3
done