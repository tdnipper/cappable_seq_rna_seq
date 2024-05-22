import os
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

basedir = os.path.dirname(os.path.realpath(__file__))
ribodepleted_dir = f"{basedir}/bbduk_reads"
out_dir = f"{basedir}/bbsplit_reads"
ref_file = f"{basedir}/myco_genome.fasta"
genome = f"{basedir}/genome/hybrid_genome.fasta"
index_dir = f"{basedir}/bbsplit_reads/index"
CPU = os.cpu_count()

files = FileHandler.get_files(ribodepleted_dir, "_R1", "_R2", file_filter = lambda x: x.endswith(".fastq.gz")and "nonrRNA" in x)

out_dir_handler = FileHandler(out_dir)
out_dir_handler.make_dir()

for key in files:
    command = f"bbsplit.sh \
                in={files[key][0]} \
                in2={files[key][1]} \
                outu1={out_dir}/{key}_R1_unmapped.fastq.gz \
                outu2={out_dir}/{key}_R2_unmapped.fastq.gz \
                ref_myco_genome={ref_file} \
                ref_hybrid_genome={genome} \
                refstats={out_dir}/{key}_refstats.txt \
                threads={CPU} \
                basename={out_dir}/{key}_#_%.fastq.gz \
                -Xmx44g"
    runner = ShellProcessRunner(command)