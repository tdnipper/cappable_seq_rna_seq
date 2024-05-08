import os
import shutil
import subprocess

basedir=os.path.abspath(os.path.dirname(__file__))
PUID = os.getuid()
PGID = os.getgid()
reads_dir=f"{basedir}/ribodepleted_reads"

if not os.path.exists(reads_dir):
    raise FileNotFoundError(f"Error: ribodepleted_reads directory does not exist at {reads_dir}")

star_dir = f"{basedir}/star_alignment"
if not os.path.exists(star_dir):
    os.mkdir(star_dir)
    os.chown(star_dir, PUID, PGID)

star_index_dir = f"{basedir}/genome/star_index"
if not os.path.exists(star_index_dir):
    raise FileNotFoundError(f"Error: STAR index directory does not exist at {star_index_dir}")

reads = {}

for dirpath, dirnames, filenames in os.walk(reads_dir):
    for filename in filenames:
        name = filename.split("_R")[0]
        if name not in reads:
            reads[name] = []
        reads[name].append(filename)

for key in reads:
    subprocess.run(["STAR",
                    "-runThreadN",
                    "4",
                    "--genomeDir",
                    star_index_dir,
                    "--readFilesIn",
                    f"{reads_dir}/{reads[key][0]}",
                    f"{reads_dir}/{reads[key][1]}",
                    "--outFileNamePrefix",
                    f"{star_dir}/{key}/{key}_",
                    "--quantMode",
                    "TranscriptomeSAM",
                    "--genomeLoad",
                    "LoadAndKeep",
                    "--outReadsUnmapped",
                    "Fastx",
                    "--readFilesCommand",
                    "zcat"])
    subprocess.run(["samtools",
                    "view",
                    f"{star_dir}/{key}/{key}_Aligned.out.sam",
                    "-o",
                    f"{star_dir}/{key}/{key}_aligned.bam"])
    os.remove(f"{star_dir}/{key}/{key}_Aligned.out.sam")
    subprocess.run(["samtools",
                    "sort",
                    f"{star_dir}/{key}/{key}_Aligned.toTranscriptome.out.bam",
                    "-o",
                    f"{star_dir}/{key}/{key}_Aligned.toTranscriptome.sorted.bam",
                    "-@",
                    "4"])
    subprocess.run(["samtools",
                    "sort",
                    f"{star_dir}/{key}/{key}_aligned.bam",
                    "-o",
                    f"{star_dir}/{key}/{key}_sorted.bam",
                    "-@",
                    "4"])
    subprocess.run(["samtools",
                    "index",
                    f"{star_dir}/{key}/{key}_Aligned.toTranscriptome.sorted.bam",
                    "-@",
                    "4"])
    subprocess.run(["samtools",
                    "index",
                    f"{star_dir}/{key}/{key}_sorted.bam",
                    "-@",
                    "4"])
    
subprocess.run(["STAR",
                "--genomeLoad Remove",
                "--genomeDir",
                star_index_dir,
                "--outFileNamePrefix",
                f"{star_dir}/exit/exit"])
shutil.rmtree(f"{star_dir}/exit")
    