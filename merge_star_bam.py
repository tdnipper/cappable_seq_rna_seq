import os
import subprocess

PUID = os.getuid()
PGID = os.getgid()
trinity_dir = "trinity_genome_guided"

if not os.path.exists(trinity_dir):
    os.makedirs(trinity_dir)
    os.chown(trinity_dir, PUID, PGID)

if not os.path.exists("star_alignment"):
    raise FileNotFoundError("Star alignment directory not at './star_alignment'")

# Get a list of bam files from star_alignment
bam_files = []
for dirpath, dirnames, filenames in os.walk("star_alignment"):
    for filename in filenames:
        if "_aligned.bam" in filename:
            bam_files.append(os.path.join(dirpath, filename))

bam_files_str = " ".join(bam_files)

# Merge bam files
result = subprocess.run(f"samtools \
                merge \
                -o {trinity_dir}/merged.bam \
                *{bam_files_str} \
                -@ 4",
                shell=True)
if result.returncode != 0:
    raise ValueError("Error merging bam files")

# Sort bam file by coordinate
subprocess.run(f"samtools \
                sort \
                -o {trinity_dir}/merged_sorted.bam \
                {trinity_dir}/merged.bam \
                -@ 4")