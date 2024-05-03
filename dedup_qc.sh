#! /bin/bash
source env.sh

sudo mkdir -p $BASEDIR/qc_reports/dedup
sudo chown $USER_ID:$GROUP_ID $BASEDIR/qc_reports/dedup
sudo chmod -R 775 $BASEDIR/qc_reports/dedup

for file in $(find $BASEDIR -name "*dedup.bam"); do
    echo "Processing $file"
    fastqc $file  -o $BASEDIR/qc_reports/dedup -t 4 --memory 10000
done