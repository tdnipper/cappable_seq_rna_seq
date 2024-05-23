import os
from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

basedir = os.path.dirname(os.path.realpath(__file__))
ribodepleted_dir = f"{basedir}/bbduk_reads"
out_dir = f"{basedir}/bbsplit_reads"
ref_file = f"{basedir}/myco_genome.fasta"
genome = f"{basedir}/genome/hybrid_genome.fasta"
index_dir = f"{basedir}/ref"
file_filter = lambda x: "nonrRNA" in x
CPU = os.cpu_count()

ribodepleted_handler = FileHandler(ribodepleted_dir)
files = ribodepleted_handler.get_files("_R1", "_R2", file_filter=file_filter)

out_dir_handler = FileHandler(out_dir)
out_dir_handler.make_dir()

index_handler = FileHandler(index_dir)
if not os.path.exists(index_dir):
    index_handler.make_dir()
    index_runner = ShellProcessRunner(f"bbsplit.sh ref_hybrid_genome={genome} ref_myco_genome={ref_file}")
    index_runner.run_shell()


for key in files:
    command = ShellProcessRunner(
                f"bbsplit.sh \
                in={files[key][0]} \
                in2={files[key][1]} \
                outu1={out_dir}/{key}_R1_unmapped.fastq.gz \
                outu2={out_dir}/{key}_R2_unmapped.fastq.gz \
                refstats={out_dir}/{key}_refstats.txt \
                threads={CPU} \
                basename={out_dir}/{key}_#_%.fastq.gz \
                -Xmx60g"
                )
    command.run_shell()
