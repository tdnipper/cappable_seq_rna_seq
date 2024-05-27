from filename_utils import FileHandler
from filename_utils import ShellProcessRunner

trimmed_dir = "bbsplit_reads"
out_dir = "nontRNA_reads"
ref_file = "genome/trna_sequences.fasta"

files_handler = FileHandler(trimmed_dir)
files = files_handler.get_files("_1", "_2", file_filter=lambda x: x.endswith("hybrid_genome.fastq.gz"))

outdir_handler = FileHandler(out_dir)
outdir_handler.make_dir()

for sample in files:
    bbduk_runner = ShellProcessRunner(f"bbduk.sh \
                                        in={files[sample][0]} \
                                        in2={files[sample][1]} \
                                        out={out_dir}/{sample}_R1_nontRNA.fastq.gz \
                                        out2={out_dir}/{sample}_R2_nontRNA.fastq.gz \
                                        outm={out_dir}/{sample}_R1_tRNA.fastq.gz \
                                        outm2={out_dir}/{sample}_R2_tRNA.fastq.gz \
                                        stats={out_dir}/{sample}_stats.txt \
                                        k=27 \
                                        ref={ref_file}")
    bbduk_runner.run_shell()