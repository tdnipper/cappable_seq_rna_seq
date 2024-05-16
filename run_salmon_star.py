#! /usr/bin/env python3

import os
import sys
import subprocess
from filename_utils import get_filenames

# Use this script to run salmon on star alignments to transcriptome
# Assumes star index has already been created

star_dir = "star_alignment"
salmon_dir = "salmon_quantification_star"
salmon_index = "trinity/salmon_index"  # Trinity generated index instead of hg38
transcripts = "genome/hybrid_transcripts_gffread.fasta"
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
        if filename.endswith("Aligned.toTranscriptome.out.bam"):
            name = filename.split("_")[0]
            print(f"Running salmon on {name}")
            result = subprocess.run(
                f"salmon quant -t {transcripts} -l A -a {dirpath}/{filename} -p 4 -o {salmon_dir}/{name}",
                shell=True,
            )
