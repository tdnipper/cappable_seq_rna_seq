#! /usr/bin/env python3

import os
import sys
import subprocess
from filename_utils import get_filenames

# Use this script to run salmon on star alignments to transcriptome
# Assumes star index has already been created

star_dir = "star_alignment"
salmon_dir = "salmon_quantification_star"
salmon_index = "trinity/salmon_index" # Trinity generated index instead of hg38
PUID = os.getuid()
PGID = os.getgid()

if not os.path.exists(salmon_dir):
    os.makedirs(salmon_dir)
    os.chown(salmon_dir, PUID, PGID)

if not os.path.exists(star_dir):
    raise FileNotFoundError("Error: star alignments directory not found")

if not os.path.exists(salmon_index):
    raise FileNotFoundError("Error: salmon index directory not found")

for dirpath, dirnames, filenames in os.walk(star_dir):
    for filename in filenames:
        if filename.endswith("Aligned.toTranscriptome.sorted.bam"):
            name = filename.strip("_Aligned.toTranscriptome.sorted.bam")
            result = subprocess.run(
                f"salmon quant -i {salmon_index} -l ISD -a {star_dir}.{filename} -p 4 -o {salmon_dir}/{name}",
                shell=True
            )