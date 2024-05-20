import os
import subprocess

leftfiles = []
rightfiles = []

reads = {}

trinity_outdir = 'trinity_independent/'
reads_dir = 'raw_data'

for dirpath, dirnames, filenames in os.walk(reads_dir):
    for filename in filenames:
        if filename.endswith("_R1.fastq.gz") or filename.endswith("_R2.fastq.gz"):
            basename = filename.rsplit("_", 1)[0]
            if basename not in reads:
                reads[basename] = {}
            if filename.endswith("_R1.fastq.gz"):
                reads[basename][0] = os.path.join(dirpath, filename)
            elif filename.endswith("_R2.fastq.gz"):
                reads[basename][1] = os.path.join(dirpath, filename)

for basename in sorted(reads.keys()):
    leftfiles.append(reads[basename][0])
    rightfiles.append(reads[basename][1])


leftfiles = ",".join(sorted(leftfiles))
rightfiles = ",".join(sorted(rightfiles))

result = subprocess.run(f"Trinity \
                        --left {leftfiles} \
                        --right {rightfiles} \
                        --seqType fq \
                        --CPU 4 \
                        --max_memory 35G \
                        --output {trinity_outdir} \
                        --trimmomatic \
                        --full_cleanup",
                        shell=True)