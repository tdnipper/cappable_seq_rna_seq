import os
from filename_utils import get_filenames_filepaths
import subprocess

basedir = os.path.dirname(os.path.realpath(__file__))
ribodepleted_dir = f"{basedir}/bbduk_reads"
out_dir = f"{basedir}/bbsplit_reads"
ref_file = f"{basedir}/myco_genome.fasta"
genome = f"{basedir}/genome/hybrid_genome.fasta"
index_dir = f"{basedir}/bbsplit_reads/index"
PUID = os.getuid()
PGID = os.getgid()
CPU = os.cpu_count()

files = get_filenames_filepaths(ribodepleted_dir, "_R1", "_R2", file_filter = lambda x: x.endswith(".fastq.gz")and "nonrRNA" in x)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    os.chown(out_dir, PUID, PGID)

for key in files:
    result = subprocess.run(f"bbsplit.sh \
                            in={files[key][0]} \
                            in2={files[key][1]} \
                            outu1={out_dir}/{key}_R1_unmapped.fastq.gz \
                            outu2={out_dir}/{key}_R2_unmapped.fastq.gz \
                            ref_myco_genome={ref_file} \
                            ref_hybrid_genome={genome} \
                            refstats={out_dir}/{key}_refstats.txt \
                            threads={CPU} \
                            basename={out_dir}/{key}_#_%.fastq.gz \
                            -Xmx44g",
                            shell=True,
                            check=True)
    if result.returncode != 0:
        raise Exception(f"Error: bbsplit failed on {key}")