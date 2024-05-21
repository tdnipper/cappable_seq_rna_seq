import os
from filename_utils import get_filenames_filepaths
import subprocess

trimmed_dir = "trimmed_reads"
out_dir = "bbduk_reads"
ref_file = "human_rRNAs.fasta"
PUID = os.getuid()
PGID = os.getgid()

files = get_filenames_filepaths(trimmed_dir, "_R1", "_R2", file_filter = lambda x: x.endswith(".fastq.gz"))

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    os.chown(out_dir, PUID, PGID)


for sample in files:
    result = subprocess.run(f"bbduk.sh \
                            in={files[sample][0]} \
                            in2={files[sample][1]} \
                            out={out_dir}/{sample}_R1_nonrRNA.fastq.gz \
                            out2={out_dir}/{sample}_R2_nonrRNA.fastq.gz \
                            outm={out_dir}/{sample}_R1_rRNA.fastq.gz \
                            outm2={out_dir}/{sample}_R2_rRNA.fastq.gz \
                            stats={out_dir}/{sample}_stats.txt \
                            k=27 \
                            ref={ref_file}",
                            shell=True,
                            check=True)
    if result.returncode != 0:
        raise Exception(f"Error: bbduk failed on {sample}")