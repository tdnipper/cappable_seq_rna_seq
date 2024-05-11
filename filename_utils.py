#! /usr/bin/env python3
__name__ = "filename_utils"
__author__ = "Thomas"
__date__ = "2024-05-08" 

import os

# Assumes prefixes are the same for read files with _R1 and _R2 denoting fwd and rev reads

def get_filenames(directory):
    replicates = {}

    for dirpath, dirnames, filenames in os.walk(directory):
        for file in filenames:
            if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
                name = file.split("_R")[0]
                if name not in replicates:
                    replicates[name] = []
                replicates[name].append(file)

    return replicates