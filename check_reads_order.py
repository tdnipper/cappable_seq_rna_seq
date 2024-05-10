import sys
import os



def check_reads_order(file1, file2):
    file1_lines = {}
    file2_lines = {}
    with open(file1, "r") as f1:
        line_number_f1 = 0
        for line in f1.readlines():
            line_number_f1 += 1
            if line.startswith("@"):
                line = line.split(" ")[0]
                file1_lines[line_number_f1] = line
    with open(file2, "r") as f2:
        line_number_f2 = 0
        for line in f2.readlines():
            line_number_f2 += 1
            if line.startswith("@"):
                line = line.split(" ")[0]
                file2_lines[line_number_f2] = line
    for key in file1_lines:
        if file1_lines[key] != file2_lines[key]:
            print(f"Mismatch detected at line {key} in {file1} and line {key} in {file2}")
            return False
    print("Files are identical")


file1 = sys.argv[1]
file2 = sys.argv[2]
check_reads_order(file1, file2)

def check_reads_order_dir(reads_dir):
    file1_lines = {}
    file2_lines = {}
    files_R1 = []
    files_R2 = []
    for dirpaths, dirnames, filenames in os.walk(reads_dir):
        for filename in filenames:
            if filename.endswith(".fastq.gz"):
                if "R1" in filename:
                    files_R1.append(filename)
                elif "R1" in filename:
                    files_R2.append(filename)
                else:
                    continue
                