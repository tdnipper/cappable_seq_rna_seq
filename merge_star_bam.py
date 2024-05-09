import os
import subprocess

PUID = os.getuid()
PGID = os.getgid()

if not os.path.exists("trinity"):
    os.makedirs("trinity")
    os.chown("trinity", PUID, PGID)

if not os.path.exists("star_alignment"):
    raise FileNotFoundError("Star alignment directory not at './star_alignment'")

# Get a list of bam files from star_alignment
bam_files = []
for dirpath, dirnames, filenames in os.walk("star_alignment"):
    for filename in filenames:
        if "_aligned.bam" in filename:
            bam_files.append(os.path.join(dirpath, filename))

# Merge bam files
subprocess.run(["samtools",
                "merge",
                "-o",
                "trinity/merged.bam",
                *bam_files
                ])

# Sort bam file by coordinate
subprocess.run(["samtools",
                "sort",
                "-o",
                "trinity/merged_sorted.bam",
                "trinity/merged.bam",
                "-@",
                "4"])