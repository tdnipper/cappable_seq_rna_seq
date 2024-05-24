from filename_utils import FileHandler
from filename_utils import ShellProcessRunner
import os

basedir = os.path.abspath(os.path.dirname(__file__))
star_dir = f"{basedir}/star_alignment"
file_filter = lambda x: "_coord_sorted.bam" in x and ".bai" not in x
outdir = f"{basedir}/qualimap/rnaseq"
gff_file = f"{basedir}/genome/hybrid_annotated_cat.gtf"
CPU = os.cpu_count()

outdir_handler = FileHandler(outdir)
outdir_handler.make_dir()

star_dir_handler = FileHandler(star_dir)
files = star_dir_handler.get_files("_coord", file_filter=file_filter)

for file in files:
    print(f"Running qualimap on {file}")
    print(f"Accessing {files[file][0]}")
    qualimap_runner = ShellProcessRunner(f"qualimap \
                                         rnaseq \
                                         -bam {files[file][0]} \
                                        -outdir {outdir}/{file} \
                                        -gtf {gff_file} \
                                        -oc {outdir}/{file}/{file}_counts.txt \
                                        -s \
                                        -pe \
                                        --java-mem-size=40G")
    qualimap_runner.run_shell()