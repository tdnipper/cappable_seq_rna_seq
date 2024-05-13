#! /usr/bin/env python3
import os
import sys
import subprocess
from filename_utils import get_filenames

__name__ = "salmon_trinity"
__author__ = "Thomas"
__date__ = "2024-05-013"

salmon_index_dir = "trinity/salmon_index"
basedir = os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
salmon_dir = f"{basedir}/salmon_quantification_trinity"
reads_dir = f"{basedir}/ribodepleted_reads/"
trinity_file = f"{basedir}/trinity/Trinity-GG.fasta"

# Use this script to run Salmon quantification on the Trinity assembled transcripts

# Check for trinity directory
if not os.path.exists("trinity"):
    raise FileNotFoundError("Error: trinity directory does not exist at ./trinity!")
else:
    os.chown("trinity", PUID, PGID)

# Check for Trinity-GG.fasta
if not os.path.exists("trinity/Trinity-GG.fasta"):
    raise FileNotFoundError("Error: Trinity-GG.fasta not found in ./trinity!")

# Check if trinity specific salmon index directory exists
# If not, create it using trinity_decoys.txt
# If trinity_decoys.txt does not exist, create it
if not os.path.exists(salmon_index_dir):
    print("trinity/salmon_index directory does not exist!")
    print("Attempting to create trinity/salmon_index directory...")
    os.mkdir(salmon_index_dir)
    os.chown(salmon_index_dir, PUID, PGID)
    # if not os.path.exists("trinity/trinity_decoys.txt"):
    #     print("trinity_decoys.txt not found, creating...")
    #     result = subprocess.run(["grep",
    #                     "'^>'",
    #                     "<(cat trinity/Trinity-GG.fasta) | cut -d ' ' -f 1 > trinity/trinity_decoys.txt"])
    #     if result.returncode != 0:
    #         raise RuntimeError("Error: Failed to create trinity_decoys.txt!")
    #     result = subprocess.run(["sed",
    #                     "-i.bak",
    #                     "-e",
    #                     "s/>//g",
    #                     "trinity/trinity_decoys.txt"])
        # if result.returncode != 0:
        #     raise RuntimeError("Error: Failed to edit trinity_decoys.txt!")
    result = subprocess.run(["salmon",
                    "index",
                    "-t",
                    trinity_file,
                    # "-d",
                    # "trinity/trinity_decoys.txt",
                    "-i",
                    salmon_index_dir,
                    "-p",
                    "4"])
    if result.returncode != 0:
        raise RuntimeError("Error: Failed to create trinity salmon index!")
    
if not os.path.exists(salmon_dir):
    os.mkdir(salmon_dir)
    os.chown(salmon_dir, PUID, PGID)

# Find all the fastq files in the ribodepleted directory
files = get_filenames(reads_dir)

# Run salmon quantification aligning to Trinity assembled transcripts
for key in files:
    if "_nonrRNA" in key:
        result = subprocess.run(["salmon",
                        "quant",
                        "-i",
                        salmon_index_dir,
                        "-l",
                        "ISF",
                        "-1",
                        f"{reads_dir}/{files[key][0]}",
                        "-2",
                        f"{reads_dir}/{files[key][1]}",
                        "-p",
                        "4",
                        "-o",
                        f"{salmon_dir}/{key}"])
        if result.returncode != 0:
            raise RuntimeError(f"Error: Failed to run salmon quantification on {key}!")
    