from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

trimmed_dir = "trimmed_reads"
out_dir = "bbduk_reads"
ref_file = "human_rRNAs.fasta"

files_handler = FileHandler(trimmed_dir)
files = files_handler.get_files()

outdir_handler = FileHandler(out_dir)
if not outdir_handler.check_exists():
    outdir_handler.make_dir()

for sample in files:
    bbduk_runner = ShellProcessRunner(f"bbduk.sh \
                                        in={files[sample][0]} \
                                        in2={files[sample][1]} \
                                        out={out_dir}/{sample}_R1_nonrRNA.fastq.gz \
                                        out2={out_dir}/{sample}_R2_nonrRNA.fastq.gz \
                                        outm={out_dir}/{sample}_R1_rRNA.fastq.gz \
                                        outm2={out_dir}/{sample}_R2_rRNA.fastq.gz \
                                        stats={out_dir}/{sample}_stats.txt \
                                        k=27 \
                                        ref={ref_file}")
    bbduk_runner.run_shell()