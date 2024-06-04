import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="Annotate gene_ids from GENCODDE with gene_name in csv file."
)
parser.add_argument(
    "files", metavar="F", type=str, nargs="+", help="File(s) to annotate."
)
args = parser.parse_args()
analyze_files = args.files


def annotate_genes(file):
    f = pd.read_csv(file)
    basename = file.split(".")[0]

    # Parse GTF file and create a dictionary mapping gene IDs to gene names
    gene_dict = {}
    with open("/home/ubuntu/blockvolume/cappable_seq_rna_seq/genome/hybrid_annotated_cat.gtf", "r") as gtf:
        for line in gtf:
            if "gene_name" in line:
                gene_id = line.split()[9].strip('"').strip('";')
                gene_name = line.split("gene_name")[1].split(";")[0].strip().strip('"')
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = gene_name
                else:
                    continue

    # Annotate genes in DataFrame
    try:
        f["gene_name"] = f.iloc[:, 0].map(gene_dict)
    except:
        f["gene_name"] = "unnamed"

    f.rename(columns={f.columns[0]: "gene_id"}, inplace=True)
    return f.to_csv(f"{basename}_annotated.csv", index=False)


for file in analyze_files:
    annotate_genes(file)
