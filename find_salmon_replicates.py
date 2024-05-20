#! /usr/bin/env python3

import os
import sys
import subprocess
from filename_utils import get_filenames_filepaths

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
reads_dir = f"{basedir}/bbduk_reads"
salmon_index_dir = f"{basedir}/genome/salmon_index"
transcript_file = f"{basedir}/genome/hybrid_transcript_gffread.fasta"
hybrid_gentrome_file = f"{salmon_dir}/hybrid_gentrome.fasta"
human_genome_file = f"{basedir}/genome/GRCh38.primary_assembly.genome.fa"
wsn_genome_file = f"{basedir}/genome/WSN_Mehle.fasta"
hybrid_genome_file = f"{basedir}/genome/hybrid_genome.fasta"

if not os.path.exists(salmon_dir):
    os.mkdir(salmon_dir)
    os.chown(salmon_dir, PUID, PGID)

if not os.path.exists(reads_dir):
    sys.exit("Error: ribodepleted_reads directory does not exist")

if not os.path.exists(salmon_index_dir):
    print("Salmon index directory does not exist")
    print("attempting salmon index")
    if not os.path.exists(f"{hybrid_gentrome_file}"):
        print("creating hybrid_gentrome.fasta")
        # Make hybrid genome file
        # if not os.path.exists(f"{basedir}/genome/hybrid_genome.fasta"):
        result = subprocess.run(
            f"cat {human_genome_file} {wsn_genome_file} > {hybrid_genome_file}",
            shell=True,
        )
        if result.returncode != 0:
            sys.exit("Error: creating hybrid_genome.fasta failed")
        # Make decoy file
        if not os.path.exists(f"{basedir}/genome/decoys.txt"):
            result = subprocess.run(
                f"cat {basedir}/genome/hybrid_genome.fasta | grep '^>' | cut -d ' ' -f 1 > {basedir}/genome/decoys.txt && \
                    sed -i.bak -e 's/>//g' {basedir}/genome/decoys.txt",
                shell=True,
            )
            if result.returncode != 0:
                sys.exit("Error: creating decoys.txt failed")
        # Make hybrid transcript file
        if not os.path.exists(f"{transcript_file}"):
            raise FileNotFoundError("hybrid_transcripts.fasta not found, please extract transcripts from hybrid genome!")
        # Make hybrid_gentrome file
        result = subprocess.run(
            f"cat {transcript_file} {hybrid_genome_file} > {hybrid_gentrome_file}",
            shell=True,
        )
        if result.returncode != 0:
            sys.exit("Error: creating hybrid_gentrome.fasta failed")
    # Run salmon index
    result = subprocess.run(
            f"salmon \
            index \
            -t {hybrid_gentrome_file} \
            -d genome/decoys.txt \
            -p 4 \
            -i {salmon_index_dir} \
            --gencode",
            shell=True
    )
    if result.returncode != 0:
        sys.exit("Error: Salmon index failed")


# Find replicate files in reads_dir and run salmon quantification on them,
# combining the replicates for each sample

# replicates = get_filenames(reads_dir)
replicates = get_filenames_filepaths(
    reads_dir, "_R1", "_R2", file_filter=lambda x: x.endswith(".fastq.gz")
)

# for key in replicates:
#     if "_nonrRNA" in key:
#         name = key.strip("_nonrRNA")
#         subprocess.run(
#             [
#                 "salmon",
#                 "quant",
#                 "-i",
#                 salmon_index_dir,
#                 "-l",
#                 "ISF",  # This might not be correct or consistent, check log files
#                 # f"<(cat {reads_dir + replicates[key][0]} {reads_dir + replicates[key][1]})"
#                 "-1",
#                 reads_dir + "/" + replicates[key][0],
#                 "-2",
#                 reads_dir + "/" + replicates[key][1],
#                 "--validateMappings",
#                 "-p",
#                 "4",
#                 "--seqBias",
#                 "--gcBias",
#                 "--reduceGCMemory",
#                 "--writeUnmappedNames",
#                 "-o",
#                 salmon_dir + "/" + name,
#             ]
#         )

for key in replicates:
    result = subprocess.run(
        f"salmon \
                            quant \
                            -i {salmon_index_dir} \
                            -l ISF \
                            -1 {replicates[key][0]} \
                            -2 {replicates[key][1]} \
                            --validateMappings \
                            -p 4 \
                            --seqBias \
                            --gcBias \
                            --reduceGCMemory \
                            --writeUnmappedNames \
                            -o {salmon_dir}/{key}",
        shell=True,
    )
    if result.returncode != 0:
        sys.exit(f"Error: salmon failed on {key}")

# Rename the output files to include the sample name

for dirpath, dirnames, filenames in os.walk(salmon_dir):
    for file in filenames:
        if "quant.sf" in file:
            sample_name = os.path.basename(dirpath)
            os.rename(
                (os.path.join(dirpath, file)),
                (os.path.join(dirpath, sample_name + "_quant.sf")),
            )
