#! /usr/bin/env python3

import os
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

# Use this script to run salmon on star alignments to transcriptome
# Assumes star index has already been created

star_dir = "star_alignment"
salmon_dir = "salmon_quantification_star"
transcripts = "genome/hybrid_exon.fasta"
PUID = os.getuid()
PGID = os.getgid()
CPU = os.cpu_count()
file_filter = lambda x: x.endswith("_Aligned.toTranscriptome.out.bam")

salmon_handler = FileHandler(salmon_dir)
salmon_handler.make_dir()
star_dir_handler = FileHandler(star_dir)
star_dir_handler.check_exists("Error: star alignments directory not found, run star first!")

files = FileHandler.get_files(star_dir_handler, "_Aligned.toTranscriptome.out.bam", file_filter=file_filter)

for name in files:
    salmon_runner = ShellProcessRunner(
        f"salmon quant -t {transcripts} -l ISR -a {files[name][0]} -p {CPU} -o {salmon_dir}/{name} --numGibbsSamples 30")
    salmon_runner.run_shell()
