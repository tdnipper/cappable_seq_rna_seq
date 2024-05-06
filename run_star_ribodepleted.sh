#!/bin/bash

source env.sh

set -e
trap cleanup EXIT

cleanup () {
    echo "STAR alignment failed, removing genome from memory"
    STAR --genomeLoad Remove \
        --genomeDir genome/star_index \
        --outFileNamePrefix star_alignments/exit/exit # remove the genome from memory
    sudo rm -R star_alignments/exit || true # remove the exit folder
    echo "STAR genome removal complete"

}

sudo mkdir -p /home/ubuntu/blockvolume/dusp11_clip-seq/ribodepleted_star_alignments
sudo chown -R $USER_ID:$GROUP_ID /home/ubuntu/blockvolume/dusp11_clip-seq/ribodepleted_star_alignments
sudo chmod -R 775 /home/ubuntu/blockvolume/dusp11_clip-seq/ribodepleted_star_alignments

# Align using TranscriptomeSAM
for file in ribodepleted_reads/*_trimmed_nonrRNA.fq.gz; do
    g=$(basename "$file" _trimmed_nonrRNA.fq.gz)
    p=${g#ultraplex_demux_}
    # d="${p}_transcriptome"
    echo "Aligning $p"
    STAR --runThreadN 4 \
        --genomeDir genome/star_index \
        --readFilesIn $file \
        --outFileNamePrefix ribodepleted_star_alignments/${p}/${p}_ribodepleted_ \
        --quantMode TranscriptomeSAM \
        --genomeLoad LoadAndKeep \
        --outReadsUnmapped Fastx \
        --readFilesCommand zcat
    echo "Converting $p to bam"
    samtools view -o ribodepleted_star_alignments/${p}/${p}_ribodepleted_aligned.bam \
        ribodepleted_star_alignments/${p}/${p}_ribodepleted_Aligned.out.sam # convert to bam
    sudo rm ribodepleted_star_alignments/${p}/${p}_ribodepleted_Aligned.out.sam # remove the sam file
    echo "Sorting $p"
    samtools sort ribodepleted_star_alignments/${p}/${p}_ribodepleted_Aligned.toTranscriptome.out.bam \
        -o ribodepleted_star_alignments/${p}/${p}_ribodepleted_transcriptome_sorted.bam \
        -@ 4 # sort the transcriptome bam file
    samtools sort ribodepleted_star_alignments/${p}/${p}_ribodepleted_aligned.bam \
        -o ribodepleted_star_alignments/${p}/${p}_ribodepleted_sorted.bam \
        -@ 4 # sort the bam file
    echo "Indexing $p"
    samtools index ribodepleted_star_alignments/${p}/${p}_ribodepleted_sorted.bam \
        -o ribodepleted_star_alignments/${p}/${p}_ribodepleted_sorted.bai \
        -@ 4 # index the sorted bam file
    samtools index ribodepleted_star_alignments/${p}/${p}_ribodepleted_transcriptome_sorted.bam \
        -o ribodepleted_star_alignments/${p}/${p}_ribodepleted_transcriptome_sorted.bai \
        -@ 4 # index the sorted transcriptome bam file
done

# Cleanup RAM
echo "STAR alignment complete, removing genome from memory"
STAR --genomeLoad Remove \
    --genomeDir genome/star_index \
    --outFileNamePrefix star_alignments/exit/exit # remove the genome from memory
sudo rm -R star_alignments/exit # remove the exit folder
echo "STAR genome removal complete"