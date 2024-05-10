import os
import subprocess
import shutil
from filename_utils import get_filenames

basedir=os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
trimmed_reads_dir = f"{basedir}/trimmed_reads" #Changed for myco removal
ribofile = f"{basedir}/decon_reads.fasta"
sort_dir = f"{basedir}/ribodepleted_reads"

if not os.path.exists(trimmed_reads_dir):
    raise FileNotFoundError(f"Error: trimmed_reads directory does not exist at {trimmed_reads_dir}")

if not os.path.exists(ribofile):
    raise FileNotFoundError(f"Error: human_rRNAs.fasta does not exist at {ribofile}")

if not os.path.exists(sort_dir):
    os.mkdir(sort_dir)
    os.chown(sort_dir, PUID, PGID)

replicates = get_filenames(trimmed_reads_dir)

for key in replicates:
    subprocess.run(["sortmerna",
                    "-ref",
                    ribofile,
                    "-reads",
                    f"{trimmed_reads_dir}/{replicates[key][0]}",
                    "-reads",
                    f"{trimmed_reads_dir}/{replicates[key][1]}",
                    "-aligned",
                    f"{sort_dir}/{key}rRNA", # label for RNA files
                    "-other",
                    f"{sort_dir}/{key}nonrRNA", # label for nonrRNA files
                    "--paired_in", # if either read aligns to rRNA, discard in aligned file
                    "--out2", # write non-rRNA neads to 2 files
                    "-fastx",
                    "--threads",
                    "4",
                    "-workdir",
                    f"{sort_dir}/tmp"])
    
    shutil.rmtree(f"{sort_dir}/tmp/kvdb")

shutil.rmtree(f"{sort_dir}/tmp")

# Rename files to match the naming convention used in the rest of the pipeline

for dirpath, dirnames, filenames in os.walk(sort_dir):
    for filename in filenames:
        if "fwd" in filename:
            new_name = filename.replace("fwd", "R1")
            shutil.move(os.path.join(dirpath, filename), os.path.join(dirpath, new_name))
        elif "rev" in filename:
            new_name = filename.replace("rev", "R2")
            shutil.move(os.path.join(dirpath, filename), os.path.join(dirpath, new_name))
        else:
            continue
