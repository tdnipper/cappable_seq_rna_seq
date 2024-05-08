import os
import subprocess
from filename_utils import get_filenames

basedir=os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
raw_dir = f"{basedir}/raw_data"
trimmed_dir = f"{basedir}/trimmed_reads"
log_dir = f"{basedir}/logs"

if not os.path.exists(f"{log_dir}/cutadapt"):
    os.makedirs(f"{log_dir}/cutadapt", mode=0o775)
    os.chown(log_dir, PUID, PGID)


if not os.path.exists(raw_dir):
    raise FileNotFoundError(f"Error: raw_data directory does not exist at {raw_dir}")

if not os.path.exists(trimmed_dir):
    os.mkdir(trimmed_dir)
    os.chown(trimmed_dir, PUID, PGID)

if not os.path.exists(f"{basedir}/fastqc/trimmed_reads"):
    os.makedirs(f"{basedir}/fastqc/trimmed_reads")
    os.chown(f"{basedir}/fastqc/trimmed_reads", PUID, PGID)

replicates = get_filenames(raw_dir)

# for dirpath, dirnames, filenames in os.walk(raw_dir):
#     for file in filenames:
#         name = file.split("_R")[0]
#         if name not in replicates:
#             replicates[name] = []
#         replicates[name].append(file)

for key in replicates:
    subprocess.run(["cutadapt",
                    f"{raw_dir}/{replicates[key][0]}",
                    f"{raw_dir}/{replicates[key][1]}",
                    "-o",
                    f"{trimmed_dir}/{key}_R1_trimmed.fastq.gz",
                    "-p",
                    f"{trimmed_dir}/{key}_R2_trimmed.fastq.gz",
                    "-q",
                    "20",
                    "--minimum-length",
                    "20",
                    "-j",
                    "3",
                    f"--json={log_dir}/cutadapt/{key}.cutadapt.json",
                    ])

subprocess.run(["fastqc",
                f"{trimmed_dir}/*",
                "-o",
                f"{basedir}/fastqc/trimmed_reads/",
                "-t",
                "3"])