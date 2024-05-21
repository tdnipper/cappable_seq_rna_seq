import os
import subprocess

leftfiles = []
rightfiles = []

CPU = os.cpu_count()

reads = {}

trinity_outdir = 'trinity_independent/'
reads_dir = 'ribodepleted_reads'

for dirpath, dirnames, filenames in os.walk(reads_dir):
    for filename in filenames:
        if filename.endswith("_nonrRNA_R1.fq.gz") or filename.endswith("_nonrRNA_R2.fq.gz"):
            basename = filename.rsplit("_", 2)[0]
            if basename not in reads:
                reads[basename] = {}
            if filename.endswith("_nonrRNA_R1.fq.gz"):
                reads[basename][0] = os.path.join(dirpath, filename)
            elif filename.endswith("_nonrRNA_R2.fq.gz"):
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
                        --CPU {CPU} \
                        --max_memory 35G \
                        --output {trinity_outdir} \
                        --full_cleanup",
                        shell=True)