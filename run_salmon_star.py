#! /usr/bin/env python3

import os
import subprocess
from filename_utils import get_filenames_filepaths
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

# Use this script to run salmon on star alignments to transcriptome
# Assumes star index has already been created

star_dir = "star_alignment"
salmon_dir = "salmon_quantification_star"
transcripts = "genome/hybrid_transcript_gffread.fasta"
PUID = os.getuid()
PGID = os.getgid()
CPU = os.cpu_count()
file_filter = lambda x: x.endswith("_Aligned.toTranscriptome.out.bam")

# if not os.path.exists(salmon_dir):
#     os.makedirs(salmon_dir)
#     os.chown(salmon_dir, PUID, PGID)

FileHandler.make_dir(salmon_dir)

FileHandler(star_dir).check_exists("Error: star alignments directory not found, run star first!")

# if not os.path.exists(star_dir):
#     raise FileNotFoundError("Error: star alignments directory not found")

# files = get_filenames_filepaths(star_dir, "_Aligned.toTranscriptome.out.bam", file_filter=file_filter)
files = FileHandler.get_files(star_dir, "Aligned.toTranscriptome.out.bam", file_filter=file_filter)

for name in files:
#     print(f"Running salmon on {name}")
#     result = subprocess.run(
#         f"salmon quant -t {transcripts} -l A -a {files[name][0]} -p {CPU} -o {salmon_dir}/{name}",
#         shell=True
#     )
#     if result.returncode != 0:
#         raise Exception(f"Error: salmon failed on {name}")
    salmon_runner = ShellProcessRunner(f"salmon quant -t {transcripts} -l A -a {files[name][0]} -p {CPU} -o {salmon_dir}/{name}")
