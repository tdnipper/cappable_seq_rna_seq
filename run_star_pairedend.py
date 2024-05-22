import os
import shutil
import signal
import sys
from filename_utils import get_filenames_filepaths
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

# This script runs star on paired end R1 R2 reads after ribodepletion

__name__ = "run_star_pairedend"
__author__ = "Thomas"
__date__ = "2024-05-07"

basedir = os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
reads_dir = f"{basedir}/bbsplit_reads"
CPU = os.cpu_count()
file_filter = lambda x: "hybrid_genome" in x
star_dir = f"{basedir}/star_alignment"
star_index_dir = f"{basedir}/genome/star_index"

# Handle CTRL+C
def signal_handler(sig, frame):
    print("\nExiting...")
    runner_remove_genome = ShellProcessRunner(f"STAR --genomeLoad Remove --genomeDir {star_index_dir} --outFileNamePrefix {star_dir}/exit/exit")
    runner_remove_genome.run_shell()
    if os.path.exists(f"{star_dir}/exit"):
        shutil.rmtree(f"{star_dir}/exit")
    print("\nGenome removed from memory")
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

# Make directories
reads_dir_handler = FileHandler(reads_dir)
reads_dir_handler.check_exists()
star_dir_handler = FileHandler(star_dir)
star_dir_handler.make_dir()
star_index_dir_handler = FileHandler(star_index_dir)
star_index_dir_handler.make_dir()

reads = reads_dir_handler.get_files(prefix1="_1", prefix2="_2", file_filter=file_filter)

for key in reads:
    if len(reads[key]) != 2:
        raise ValueError(f"Error: {key} does not have exactly 2 read files")

# Run star alignment on each sample in paired-end mode
for key in reads:
    keyname = key  # used to split key, not needed now but keeping syntax the same
    runner_star = ShellProcessRunner(f"STAR \
                                --runThreadN {CPU} \
                                --genomeDir {star_index_dir} \
                                --readFilesIn {reads[key][0]} {reads[key][1]} \
                                --outFileNamePrefix {star_dir}/{keyname}/{keyname}_ \
                                --quantMode TranscriptomeSAM \
                                --genomeLoad LoadAndKeep \
                                --outReadsUnmapped Fastx \
                                --readFilesCommand zcat")
    runner_star.run_shell()
    
    # Convert sam to bam
    runner_convert = ShellProcessRunner(f"samtools view {star_dir}/{keyname}/{keyname}_Aligned.out.sam -o {star_dir}/{keyname}/{keyname}_aligned.bam")
    runner_convert.run_shell()
    os.remove(f"{star_dir}/{keyname}/{keyname}_Aligned.out.sam")
    
    # Sort and index bam files
    runner_sort_transcriptome = ShellProcessRunner(f"samtools sort {star_dir}/{keyname}/{keyname}_Aligned.toTranscriptome.out.bam -o {star_dir}/{keyname}/{keyname}_Aligned.toTranscriptome.sorted.bam -@ {CPU}")
    runner_sort_transcriptome.run_shell()
    runner_sort_genome = ShellProcessRunner(f"samtools sort {star_dir}/{keyname}/{keyname}_aligned.bam -o {star_dir}/{keyname}/{keyname}_coord_sorted.bam -@ {CPU}")
    runner_sort_genome.run_shell()
    runner_index_transcriptome = ShellProcessRunner(f"samtools index {star_dir}/{keyname}/{keyname}_Aligned.toTranscriptome.sorted.bam -@ {CPU}")
    runner_index_transcriptome.run_shell()
    runner_index_genome = ShellProcessRunner(f"samtools index {star_dir}/{keyname}/{keyname}_coord_sorted.bam -@ {CPU}")
    runner_index_genome.run_shell()

    
# Remove genome from memory
runner_remove_genome = ShellProcessRunner(f"STAR --genomeLoad Remove --genomeDir {star_index_dir} --outFileNamePrefix {star_dir}/exit/exit")
runner_remove_genome.run_shell()

# Clean up
if os.path.exists(f"{star_dir}/exit"):
    shutil.rmtree(f"{star_dir}/exit")
