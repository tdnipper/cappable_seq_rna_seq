import os
import subprocess
import shutil

basedir=os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
trimmed_reads_dir = f"{basedir}/trimmed_reads/"
ribofile = f"{basedir}/human_rRNAs.fasta"
sort_dir = f"{basedir}/ribodepleted_reads"

if not os.path.exists(trimmed_reads_dir):
    raise FileNotFoundError(f"Error: trimmed_reads directory does not exist at {trimmed_reads_dir}")

if not os.path.exists(ribofile):
    raise FileNotFoundError(f"Error: human_rRNAs.fasta does not exist at {ribofile}")

if not os.path.exists(sort_dir):
    os.mkdir(sort_dir)
    os.chown(sort_dir, PUID, PGID)

replicates = {}

for dirpath, dirnames, filenames in os.walk(trimmed_reads_dir):
    for file in filenames:
        name = file.split("R")[0]
        if name not in replicates:
            replicates[name] = []
        replicates[name].append(file)

for key in replicates:
    # key = key.strip("_")
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
                    "ribodepleted_reads/tmp"])
    
    shutil.rmtree("ribodepleted_reads/tmp/kvdb")

shutil.rmtree("ribodepleted_reads/tmp")

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
