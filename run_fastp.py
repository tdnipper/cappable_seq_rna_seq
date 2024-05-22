from filename_utils import FileHandler
from filename_utils import ShellProcessRunner
import os

input_dir = "raw_data"
output_dir = "trimmed_reads"
fastqc_dir = "fastqc/trimmed_reads"
CPU = os.cpu_count()

output_dir_handler = FileHandler(output_dir)
output_dir_handler.make_dir()

input_handler = FileHandler(input_dir)
files = FileHandler.get_files(input_handler, "_R1", "_R2")

for key in files:
    fastp_runner = ShellProcessRunner(f"fastp -i {files[key][0]} \
        -I {files[key][1]} \
        -o {output_dir}/{key}_R1_trimmed.fastq.gz \
        -O {output_dir}/{key}_R2_trimmed.fastq.gz \
        -j {output_dir}/{key}.fastp.json \
        -h {output_dir}/{key}.fastp.html \
        --unpaired1 {output_dir}/{key}_R1_unpaired.fastq.gz \
        --unpaired2 {output_dir}/{key}_R2_unpaired.fastq.gz \
        -l 25 \
        -q 20")
    fastp_runner.run_shell()

fastqc_runner = ShellProcessRunner(f"fastqc {output_dir}/ -o {fastqc_dir}/ -t {CPU}")
fastqc_runner.run_shell()