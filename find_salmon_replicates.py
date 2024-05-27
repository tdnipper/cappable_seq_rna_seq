#! /usr/bin/env python3

import os
import sys
import subprocess
from filename_utils import get_filenames_filepaths
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

__name__ = "find_salmon_replicates"
__author__ = "Thomas"
__date__ = "2024-05-07"

# Use this file to find all replicates denoted by _R(rep number)_ for each sample
# and run Salmon quantification on them. This script will output a directory
# containing the Salmon output files for each sample with replicates combined.

# Define the path to the directory containing the Salmon output files
basedir = os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
salmon_dir = f"{basedir}/salmon_quantification"
reads_dir = f"{basedir}/bbsplit_reads"
salmon_index_dir = f"{basedir}/genome/salmon_index"
transcript_file = f"{basedir}/genome/hybrid_exon.fasta"
hybrid_gentrome_file = f"{salmon_dir}/hybrid_gentrome.fasta"
human_genome_file = f"{basedir}/genome/GRCh38.primary_assembly.genome.fa"
wsn_genome_file = f"{basedir}/genome/WSN_Mehle.fasta"
hybrid_genome_file = f"{basedir}/genome/hybrid_genome.fasta"
CPU = os.cpu_count()

# if not os.path.exists(salmon_dir):
#     print("Creating salmon_quantification directory")
#     os.mkdir(salmon_dir)
#     os.chown(salmon_dir, PUID, PGID)

salmon_dir_handler = FileHandler(salmon_dir)
salmon_dir_handler.create_dir()

# if not os.path.exists(reads_dir):
#     sys.exit("Error: ribodepleted_reads directory does not exist")

reads_dir_handler = FileHandler(reads_dir)
reads_dir_handler.check_exists("Error: ribodepleted_reads directory does not exist")

if not os.path.exists(salmon_index_dir):
    print("Salmon index directory does not exist")
    print("attempting salmon index")
    if not os.path.exists(f"{hybrid_gentrome_file}"):
        print("creating hybrid_gentrome.fasta")
        # Make hybrid genome file
        hybrid_genome_runner = ShellProcessRunner(f"cat {human_genome_file} {wsn_genome_file} > {hybrid_genome_file}")
        hybrid_genome_runner.run_shell()

        # Make decoy file
        if not os.path.exists(f"{basedir}/genome/decoys.txt"):
            decoys_runner = ShellProcessRunner(f"cat {basedir}/genome/hybrid_genome.fasta | grep '^>' | cut -d ' ' -f 1 > {basedir}/genome/decoys.txt && \
                sed -i.bak -e 's/>//g' {basedir}/genome/decoys.txt")
            decoys_runner.run_shell()

        # Check for hybrid transcript file
        if not os.path.exists(f"{transcript_file}"):
            raise FileNotFoundError("hybrid_transcripts.fasta not found, please extract transcripts from hybrid genome!")
        
        # Make hybrid_gentrome file
        hybrid_gentrome_runner = ShellProcessRunner(f"cat {transcript_file} {hybrid_genome_file} > {hybrid_gentrome_file}")
        hybrid_gentrome_runner.run_shell()

    # Run salmon index
    salmon_index_runner = ShellProcessRunner(f"salmon index -t {hybrid_gentrome_file} -d genome/decoys.txt -p {CPU} -i {salmon_index_dir} --gencode -k 25")
    salmon_index_runner.run_shell()


# Find replicate files in reads_dir and run salmon quantification on them,
# combining the replicates for each sample
replicates_handler = FileHandler(reads_dir)
replicates = replicates_handler.get_files("_1", "_2", file_filter=lambda x: "hybrid_genome" in x)

for key in replicates:
    salmon_runner = ShellProcessRunner(f"salmon \
                        quant \
                        -i {salmon_index_dir} \
                        -l ISR \
                        -1 {replicates[key][0]} \
                        -2 {replicates[key][1]} \
                        --validateMappings \
                        -p {CPU} \
                        --seqBias \
                        --gcBias \
                        --reduceGCMemory \
                        --writeUnmappedNames \
                        -o {salmon_dir}/{key} \
                        --recoverOrphans \
                        --numGibbsSamples 30")
    salmon_runner.run_shell()

# Rename the output files to include the sample name

for dirpath, dirnames, filenames in os.walk(salmon_dir):
    for file in filenames:
        if "quant.sf" in file:
            sample_name = os.path.basename(dirpath)
            os.rename(
                (os.path.join(dirpath, file)),
                (os.path.join(dirpath, sample_name + "_quant.sf")),
            )
            