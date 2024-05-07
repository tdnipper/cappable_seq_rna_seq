#! /usr/bin/env python3

import os
import sys
import subprocess

__name__ = "find_salmon_replicates"
__author__ = "Thomas"
__date__ = "2024-05-07"

# Use this file to find all replicates denoted by _R(rep number)_ for each sample
# and run Salmon quantification on them. This script will output a directory
# containing the Salmon output files for each sample with replicates combined.

# Define the path to the directory containing the Salmon output files
salmon_dir = "salmon_quantification"
reads_dir = "ribodepleted_reads/"
salmon_index_dir = "genome/salmon_index"

PUID = os.getuid()
PGID = os.getgid()

if not os.path.exists(salmon_dir):
    os.mkdir(salmon_dir)
    os.chown(salmon_dir, PUID, PGID)

if not os.path.exists(reads_dir):
    sys.exit("Error: ribodepleted_reads directory does not exist")

if not os.path.exists(salmon_index_dir):
    print("Error: Salmon index directory does not exist")
    print("attempting salmon index")
    subprocess.run(["salmon", 
                    "index", 
                    "-t",
                    "genome/hybrid_gentrome.fasta",
                    "-d", 
                    "genome/decoys.txt",
                    "-p",
                    "4", 
                    "-i", 
                    salmon_index_dir])


# Find replicate files in reads_dir and run salmon quantification on them, 
# combining the replicates for each sample

replicates = {}
for dirpath, dirnames, filenames in os.walk(reads_dir):
    for file in filenames:
        if "nonrRNA" in file:
            name = file.split("R")[0]
            if name not in replicates:
                replicates[name] = []
            replicates[name].append(file)

for key in replicates:
    subprocess.run(["salmon",
                    "quant",
                    "-i",
                    salmon_index_dir,
                    "-l",
                    "A",
                    "-r",
                    reads_dir + replicates[key][0] + " " + reads_dir + replicates[key][1],
                    "--validateMappings",
                    "-p",
                    "4",
                    "--seqBias",
                    "--gcBias",
                    "--reduceGCMemory",
                    "--writeUnmappedNames",
                    "-o",
                    salmon_dir + "/" + key.strip("_")])
    
# Rename the output files to include the sample name

for dirpath, dirnames, filenames in os.walk(salmon_dir):
    for file in filenames:
        if "quant.sf" in file:
            sample_name = os.path.basename(dirpath)
            os.rename((os.path.join(dirpath, file)),(os.path.join(dirpath, sample_name + "_quant.sf")))