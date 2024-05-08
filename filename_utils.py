__name__ = "filename_utils"
__author__ = "Thomas"
__date__ = "2024-05-08" 

import os

def get_filenames(directory):
    replicates = {}

    for dirpath, dirnames, filenames in os.walk(directory):
        for file in filenames:
            name = file.split("_R")[0]
            if name not in replicates:
                replicates[name] = []
            replicates[name].append(file)

    return replicates