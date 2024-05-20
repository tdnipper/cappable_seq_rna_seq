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
            if ".log" not in file:
                name = file.split("_R")[0]
                if name not in replicates:
                    replicates[name] = {}
                if "R1" in file:
                    replicates[name][0] = os.path.join(dirpath, file)
                elif "R2" in file:
                    replicates[name][1] = os.path.join(dirpath, file)

    return replicates