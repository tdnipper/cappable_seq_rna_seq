#! /bin/bash

source env.sh

sudo mkdir qc_reports/qualimap
sudo chown -R $USER_ID:$GROUP_ID qc_reports/qualimap
sudo chmod -R 775 qc_reports/qualimap

function handle_interrupt() {
    echo " Cancelled"
    exit 1
}
trap handle_interrupt SIGINT

qualimap multi-bamqc -d qualimap_input_files.tsv \
    -gff genome/genomic.gtf \
    -c \
    -outdir qc_reports/qualimap/ \
    -outfile qualimap_report.html \
    -outformat HTML \
    -r \
    --java-mem-size=4G