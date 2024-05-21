import os
import shutil
import subprocess
from filename_utils import get_filenames_filepaths

# This script runs star on paired end R1 R2 reads after ribodepletion

__name__ = "run_star_pairedend"
__author__ = "Thomas"
__date__ = "2024-05-07"

basedir = os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
reads_dir = f"{basedir}/bbduk_reads"
CPU = os.cpu_count()

if not os.path.exists(reads_dir):
    raise FileNotFoundError(
        f"Error: ribodepleted_reads directory does not exist at {reads_dir}"
    )

star_dir = f"{basedir}/star_alignment"
if not os.path.exists(star_dir):
    os.mkdir(star_dir)
    os.chown(star_dir, PUID, PGID)

star_index_dir = f"{basedir}/genome/star_index"
if not os.path.exists(star_index_dir):
    raise FileNotFoundError(
        f"Error: STAR index directory does not exist at {star_index_dir}"
    )

# reads = get_filenames(reads_dir)
reads = get_filenames_filepaths(
    reads_dir, "_R1", "_R2", file_filter=lambda x: "nonrRNA" in x
)


for key in reads:
    if len(reads[key]) != 2:
        raise ValueError(f"Error: {key} does not have exactly 2 read files")

# Run star alignment on each sample in paired-end mode
for key in reads:
    keyname = key  # used to split key, not needed now but keeping syntax the same
    result = subprocess.run(
        f"STAR \
                --runThreadN {CPU} \
                --genomeDir {star_index_dir} \
                --readFilesIn {reads[key][0]} {reads[key][1]} \
                --outFileNamePrefix {star_dir}/{keyname}/{keyname}_ \
                --quantMode TranscriptomeSAM \
                --genomeLoad LoadAndKeep \
                --outReadsUnmapped Fastx \
                --readFilesCommand zcat",
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Error: STAR failed to align {key}")
    # Sort and index bam files
    result = subprocess.run(
        f"samtools view {star_dir}/{keyname}/{keyname}_Aligned.out.sam -o {star_dir}/{keyname}/{keyname}_aligned.bam",
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Error: samtools failed to convert {key} to bam")
    os.remove(f"{star_dir}/{keyname}/{keyname}_Aligned.out.sam")
    result = subprocess.run(
        f"samtools sort {star_dir}/{keyname}/{keyname}_Aligned.toTranscriptome.out.bam -o {star_dir}/{keyname}/{keyname}_Aligned.toTranscriptome.sorted.bam -@ {CPU}",
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Error: samtools failed to sort {key}_Aligned.toTranscriptome.out.bam"
        )
    result = subprocess.run(
        f"samtools sort {star_dir}/{keyname}/{keyname}_aligned.bam -o {star_dir}/{keyname}/{keyname}_coord_sorted.bam -@ {CPU}",
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Error: samtools failed to sort {key}_aligned.bam")
    result = subprocess.run(
        f"samtools index {star_dir}/{keyname}/{keyname}_Aligned.toTranscriptome.sorted.bam -@ {CPU}",
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Error: samtools failed to index {key}_Aligned.toTranscriptome.sorted.bam"
        )
    result = subprocess.run(
        f"samtools index {star_dir}/{keyname}/{keyname}_coord_sorted.bam -@ {CPU}",
        shell=True,
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Error: samtools failed to index {key}_coord_sorted.bam")

# Remove genome from memory
result = subprocess.run(
    f"STAR --genomeLoad Remove --genomeDir {star_index_dir} --outFileNamePrefix {star_dir}/exit/exit",
    shell=True,
    check=True,
)
if result.returncode != 0:
    raise RuntimeError(f"Error: STAR failed to remove genome from memory")

# Clean up
if os.path.exists(f"{star_dir}/exit"):
    shutil.rmtree(f"{star_dir}/exit")
