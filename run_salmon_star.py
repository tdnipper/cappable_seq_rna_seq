#! /usr/bin/env python3

import os
import sys
import subprocess
from filename_utils import get_filenames

# Use this script to run salmon on star alignments to transcriptome

star_dir = "star_alignments"
salmon_dir = "salmon_quantification_star"
salmon_index = "trinity/salmon_index" # Trinity generated index instead of hg38
PUID = os.getuid()
PGID = os.getgid()

if not os.path.exists(salmon_dir):
    os.makedirs(salmon_dir)
    os.chown(salmon_dir, PUID, PGID)

if not os.path.exists(star_dir):
    raise FileNotFoundError("Error: star alignments directory not found")

files = get_filenames(star_dir)

for key in files:
    if "_nonrRNA" in key:
        name=key.strip("_nonrRNA")
        result = subprocess.run(
            [
                "salmon",
                "quant",
                "-i",
                salmon_index,
                "-l",
                "ISD", #Maybe change this?
                "-1",
                f"{star_dir}/{files[key][0]}",
                "-2",
                f"{star_dir}/{files[key][1]}",
                "-p",
                "8",
                "-o",
                f"{salmon_dir}/{name}",
            ]
        )