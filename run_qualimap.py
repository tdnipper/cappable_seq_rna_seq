from filename_utils import FileHandler
from filename_utils import ShellProcessRunner
import os

basedir = os.path.abspath(os.path.dirname(__file__))
star_dir = f"{basedir}/star_alignment"
file_filter = lambda x: "_coord_sorted.bam" in x
outdir = f"{basedir}/qualimap/bamqc"
gff_file = f"{basedir}/genome/hybrid_annotated_agat.gtf"
CPU = os.cpu_count()

outdir_handler = FileHandler(outdir)
outdir_handler.make_dir()

star_dir_handler = FileHandler(star_dir)
files = star_dir_handler.get_files("_coord", file_filter=file_filter)

for file in files:
    print(f"Running qualimap on {file}")
    print(f"Accessing {files[file][0]}")
    qualimap_runner = ShellProcessRunner(f"qualimap \
                                         bamqc \
                                         -bam {files[file][0]} \
                                        -outdir {outdir}/{file} \
                                        -gff {gff_file} \
                                        -nt {CPU} \
                                        -oc {outdir}/{file}/{file}_genome_coverage.txt \
                                        -os {outdir}/{file}/{file}_stats.txt \
                                        -c \
                                        --java-mem-size=40G")
    qualimap_runner.run_shell()