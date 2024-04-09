#!/usr/bin/env python3
__author__ = "Thomas Nipper"
__version__ = "0.1.0"

import argparse

# Genbank to GTF conversion script
# GTF needs: seqname, source, feature, start, end, score, strand, frame, attribute


args = argparse.ArgumentParser(
    description="Convert Genbank to GTF format giving a list of Genbank files and a GTF file prefix"
)
args.add_argument("-g", "--genbank_file", help="Genbank file to convert", nargs="+")
args.add_argument(
    "-t", "--transcript_file", help="Transcript file to convert", nargs="+"
)
args.add_argument("-o", "--gtf_file", help="GTF file prefix to write to")
args = args.parse_args()
gtf_file = args.gtf_file
gtf_file = f"{gtf_file}.gtf"

with open(gtf_file, "w") as gtf:
    for genbank_file in args.genbank_file:
        with open(genbank_file) as f:
            lines = f.readlines()

        for line in lines:
            if "/gene" in line:
                seqname = line.split("=")[1].strip().replace('"', "")
                gene_id = f"WSN_{seqname}"
                break

        for line in lines:
            if "VERSION" in line:
                source = line.split()[1].strip()
                break

        for line in lines:
            if "source" in line:
                start, end = line.split()[1].strip().split("..")
                break
        
        for line in lines:
            if "product" in line:
                gene_type = line.split("=")[1].strip().replace('"', "")
                break
        # Write info to GTF file once per GenBank file
        gtf.write(
            f'{gene_id}\t{source}\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{gene_id}";\tgene_name "{seqname}";\tgene_biotype "{gene_type}";\n'
        )

# Append transcript info to GTF file as exons
# NOTE: I had to custom format the >line in the transcript file to get the necessary info
# This is not a general solution, but it works for the given transcript file because we have so few genes
        
with open(gtf_file, "a") as gtf:
    for transcript_file in args.transcript_file:
        with open(transcript_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if ">" in line:
                    seqname = line.split("|")[1].split (" ")[0].strip()
                    gene_id = f"WSN_{seqname}"
                    source = line.split(":")[0].strip()
                    source = source.replace(">", "")
                    start = line.split(":")[1].split("..")[0].strip()
                    end = line.split(":")[1].split("..")[1].split(" ")[0].strip()
                    gene_type = line.split("=")[1].strip()
                    # Write info to GTF file once per transcript file line that begins with ">"
                    gtf.write(
                        f'{gene_id}\t{source}\texon\t{start}\t{end}\t.\t+\t1\tgene_id "{gene_id}";\tgene_name "{seqname}";\ttranscript_id "{gene_id}_{seqname}";\ttranscript_name "{gene_type}";\tgene_biotype "{gene_type}";\n'
                    )
