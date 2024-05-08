#! /bin/bash

source env.sh

sudo mkdir ribodepleted_reads/
sudo chown $USER_ID:$GROUP_ID ribodepleted_reads
sudo chmod -R 775 ribodepleted_reads

# Sort rRNA reads
for file in $(find trimmed_reads/ -name "*_R1_trimmed.fastq.gz"); do
    base=$(basename $file ".fastq.gz")
    echo "Sorting $base"
    sortmerna -ref human_rRNAs.fasta -reads $file --threads 4 -fastx -workdir ribodepleted_reads/$base -aligned ribodepleted_reads/$base"_rRNA" -other ribodepleted_reads/$base"_nonrRNA"
    sudo rm -R ribodepleted_reads/$base/ # Cleanup working directory
done
echo "Sortmerna complete"