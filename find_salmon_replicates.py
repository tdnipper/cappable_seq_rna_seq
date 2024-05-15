#! /usr/bin/env python3

import os
import sys
import subprocess
from filename_utils import get_filenames

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
reads_dir = f"{basedir}/ribodepleted_reads"
salmon_index_dir = f"{basedir}/genome/salmon_index"

if not os.path.exists(salmon_dir):
    os.mkdir(salmon_dir)
    os.chown(salmon_dir, PUID, PGID)

if not os.path.exists(reads_dir):
    sys.exit("Error: ribodepleted_reads directory does not exist")

if not os.path.exists(salmon_index_dir):
    print("Salmon index directory does not exist")
    print("attempting salmon index")
    if not os.path.exists(f"{basedir}/genome/hybrid_gentrome.fasta"):
        print("creating hybrid_gentrome.fasta")
        # Make hybrid genome file
        if not os.path.exists(f"{basedir}/genome/hybrid_genome.fasta"):
            result = subprocess.run(
                f"cat {basedir}/genome/GRCh38.p14.genome.fa {basedir}/genome/WSN_Mehle.fasta > {basedir}/genome/hybrid_genome.fasta",
                shell=True
            )
            if result.returncode != 0:
                sys.exit("Error: creating hybrid_genome.fasta failed")
        # Make decoy file
        if not os.path.exists(f"{basedir}/genome/decoys.txt"):
            result = subprocess.run(
                f"cat {basedir}/genome/hybrid_genome.fasta | grep '^>' | cut -d ' ' -f 1 > {basedir}/genome/decoys.txt && \
                    sed -i.bak -e 's/>//g' {basedir}/genome/decoys.txt",
                shell=True
            )
            if result.returncode != 0:
                sys.exit("Error: creating decoys.txt failed")
        # Make hybrid transcript file
        if not os.path.exists(f"{basedir}/genome/hybrid_transcripts.fasta"):
            result = subprocess.run(
                f"cat {basedir}/genome/gencode.v46.transcripts.fa {basedir}/genome/WSN_transcripts.fasta > {basedir}/genome/hybrid_transcripts.fasta",
                shell=True
            )
            if result.returncode != 0:
                sys.exit("Error: creating hybrid_transcripts.fasta failed")
        # Make hybrid_gentrome file
        result = subprocess.run(
            f"cat {basedir}/genome/hybrid_genome.fasta {basedir}/genome/hybrid_transcripts.fasta > {basedir}/genome/hybrid_gentrome.fasta",
            shell=True
        )
    # Run salmon index
    result = subprocess.run(
        [
            "salmon",
            "index",
            "-t",
            "genome/hybrid_gentrome.fasta",
            "-d",
            "genome/decoys.txt",
            "-p",
            "4",
            "-i",
            salmon_index_dir,
        ]
    )
    if result.returncode != 0:
        sys.exit("Error: Salmon index failed")


# Find replicate files in reads_dir and run salmon quantification on them,
# combining the replicates for each sample

replicates = get_filenames(reads_dir)

for key in replicates:
    if "_nonrRNA" in key:
        name = key.strip("_nonrRNA")
        subprocess.run(
            [
                "salmon",
                "quant",
                "-i",
                salmon_index_dir,
                "-l",
                "ISF",  # This might not be correct or consistent, check log files
                # f"<(cat {reads_dir + replicates[key][0]} {reads_dir + replicates[key][1]})"
                "-1",
                reads_dir + "/" + replicates[key][0],
                "-2",
                reads_dir + "/" + replicates[key][1],
                "--validateMappings",
                "-p",
                "4",
                "--seqBias",
                "--gcBias",
                "--reduceGCMemory",
                "--writeUnmappedNames",
                "-o",
                salmon_dir + "/" + name,
            ]
        )

# Rename the output files to include the sample name

for dirpath, dirnames, filenames in os.walk(salmon_dir):
    for file in filenames:
        if "quant.sf" in file:
            sample_name = os.path.basename(dirpath)
            os.rename(
                (os.path.join(dirpath, file)),
                (os.path.join(dirpath, sample_name + "_quant.sf")),
            )
