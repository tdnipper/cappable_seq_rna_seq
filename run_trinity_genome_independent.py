import os
from filename_utils import ShellProcessRunner
from filename_utils import FileHandler

CPU = os.cpu_count()
trinity_outdir = 'trinity_independent/'
reads_dir = 'bbsplit_reads'

files_handler = FileHandler(reads_dir)
files = files_handler.get_files('_1', '_2', file_filter = lambda x: "hybrid_genome" in x)
left_files = ",".join(files[key][0] for key in sorted(files.keys()))
right_files = ",".join(files[key][1] for key in sorted(files.keys()))

trinity_runner = ShellProcessRunner(f"Trinity \
                        --left {left_files} \
                        --right {right_files} \
                        --seqType fq \
                        --CPU {CPU} \
                        --max_memory 35G \
                        --output {trinity_outdir} \
                        --full_cleanup \
                        --SS_lib_type RF")
trinity_runner.run_shell()
