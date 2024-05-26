#! /usr/bin/env python3
__name__ = "filename_utils"
__author__ = "Thomas"
__date__ = "2024-05-08" 

import os
import subprocess

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
                    replicates[name][0] = file
                elif "R2" in file:
                    replicates[name][1] = file

    return replicates

def get_filenames_filepaths(directory: str, prefix1: str = "_R1", prefix2: str = "_R2", file_filter: callable = None) -> dict:
    """
    Retrieves the filenames and filepaths from a given directory.

    Args:
        directory (str): The directory to search for files.
        prefix1 (str, optional): The prefix for the first set of files. Defaults to "_R1", but can be changed to 'fwd' for example.
        prefix2 (str, optional): The prefix for the second set of files. Defaults to "_R2", but can be changed to 'rev' for example.
        file_filter (function, optional): A function to filter files. Defaults to None. Can be used to filter out log files or non fastq, etc.

    Returns:
        dict: A dictionary containing the sample names and their corresponding filepaths.
              The keys of the dictionary are the filenames split at the prefixes, and the values are nested dictionaries.
              The nested dictionaries have keys 0 and 1, representing the first and second set of files respectively.
              The values of the nested dictionaries are the filepaths.
    
    Example:
        To find filenames with prefixes R1 and R2 in the directory 'trimmed_reads', excluding log files:
        get_filenames_filepaths('trimmed_reads', '_R1', '_R2', file_filter = lambda x: '.log' not in x)

    """
    files = {}
    for dirpath, dirnames, filenames in os.walk(directory):
        for file in filenames:
            if file_filter is None or file_filter(file):
                if prefix1 in file:
                    name = file.split(prefix1)[0]
                elif prefix2 in file:
                    name = file.split(prefix2)[0]
                if name not in files:
                    files[name] = {}
                if prefix1 in file:
                    files[name][0] = os.path.join(dirpath, file)
                elif prefix2 in file:
                    files[name][1] = os.path.join(dirpath, file)
                
    return files

class FileHandler:
    def __init__(self, path):
        self.path = path
    
    def check_exists(self, error_message: str = None):
        if not os.path.exists(self.path):
            if error_message is not None:
                raise FileNotFoundError(error_message)
            else:
                raise FileNotFoundError(f"Error: {self.path} does not exist.")
        return os.path.exists(self.path)
    
    def make_dir(self):
        try:
            os.makedirs(self.path, exist_ok=True)
            os.chown(self.path, os.getuid(), os.getgid())
        except Exception as e:
            raise RuntimeError(f"Error: Could not create directory {self.path}. {str(e)}")
    
    def get_files(self, prefix1: str = "_R1", prefix2: str = "_R2", file_filter: callable = None):
        """
        Retrieves the filenames and filepaths from a given directory.

        Args:
            directory (str): The directory to search for files.
            prefix1 (str, optional): The prefix for the first set of files. Defaults to "_R1", but can be changed to 'fwd' for example.
            prefix2 (str, optional): The prefix for the second set of files. Defaults to "_R2", but can be changed to 'rev' for example.
            file_filter (function, optional): A function to filter files. Defaults to None. Can be used to filter out log files or non fastq, etc.

        Returns:
            dict: A dictionary containing the sample names and their corresponding filepaths.
                The keys of the dictionary are the filenames split at the prefixes, and the values are nested dictionaries.
                The nested dictionaries have keys 0 and 1, representing the first and second set of files respectively.
                The values of the nested dictionaries are the filepaths.
        
        Example:
            To find filenames with prefixes R1 and R2 in the directory 'trimmed_reads', excluding log files:
            get_filenames_filepaths('trimmed_reads', '_R1', '_R2', file_filter = lambda x: '.log' not in x)

        """
        files = {}
        for dirpath, _, filenames in os.walk(self.path):
            for file in filenames:
                name = None
                if file_filter is None or file_filter(file):
                    if prefix1 in file:
                        name = file.split(prefix1)[0]
                    elif prefix2 in file:
                        name = file.split(prefix2)[0]
                    if name not in files:
                        files[name] = {}
                    if prefix1 in file:
                        files[name][0] = os.path.join(dirpath, file)
                    elif prefix2 in file:
                        files[name][1] = os.path.join(dirpath, file)
                
        return files

class ShellProcessRunner:
    def __init__(self, command):
        self.command = command

    def run_shell(self):
        try:
            subprocess.run(self.command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error: {e}")