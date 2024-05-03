#! /bin/zsh
sudo mkdir trimmed_reads/
sudo chown ubuntu:ubuntu trimmed_reads/
sudo chmod 775 trimmed_reads/

mkdir qc_reports/cutadapt/
sudo chown ubuntu:ubuntu qc_reports/cutadapt/
sudo chmod 775 qc_reports/cutadapt/

for i in raw_data/*.fastq.gz; do
    p=$(basename "$i" .fastq.gz)
    echo "Processing ${p}"
    cutadapt -o "trimmed_reads/${p%}_trimmed.fastq.gz" $i -q 20 --minimum-length 20 -n 3 -j 0 1> "qc_reports/cutadapt/${p%}_cutadapt.txt"
done