# This script creates a star index from the hybrid genome using hybrid gtf from agat created hybrid gff

basedir=/home/ubuntu/blockvolume/cappable_seq_rna_seq
genome_dir=genome
index_dir=genome/star_index
genome_file=genome/hybrid_genome.fasta
gff_file=genome/hybrid_annotated_agat_sort.gtf
USER_ID=$(id -u)
GROUP_ID=$(id -g)

# Set memory limit (in kilobytes, here it is 38GB)
ulimit -v 38000000

if [ ! -d $basedir/$index_dir ]; then
    mkdir -p $basedir/$index_dir
    chown $USER_ID:$GROUP_ID $basedir/$index_dir
fi

if [ ! -f $basedir/$genome_file ]; then
    echo "Genome file not found: $basedir/$genome_file"
    exit 1
fi

STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir $basedir/$index_dir \
    --genomeFastaFiles $basedir/$genome_file \
    --sjdbGTFfile $basedir/$gff_file \
    --sjdbOverhang 149 \
    --limitGenomeGenerateRAM 36000000000 \
    --sjdbGTFtagExonParentTranscript Parent \

# For 150bp paired end sequencing, the overhang should be read length - 1


