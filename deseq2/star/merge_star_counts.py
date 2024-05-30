import pandas as pd
import os

star_dir = "/home/ubuntu/blockvolume/cappable_seq_rna_seq/star_alignment"

def combine_star_counts(directory):
    data = pd.DataFrame()
    for dirpath, dirnames, filenames in os.walk(directory):
        for file in filenames:
            if file.endswith("_ReadsPerGene.out.tab"):
                basename = file.split("_")[0]
                raw_data = pd.read_csv(os.path.join(dirpath, file), sep="\t")
                colnames = raw_data.columns.tolist()
                colnames[0] = "gene_id"
                colnames[1] = "unstranded"
                colnames[2] = "forward"
                colnames[3] = "reverse"
                raw_data.columns = colnames
                raw_data = raw_data[["gene_id", "reverse"]].iloc[4:]
                if data.empty:
                    data = pd.concat([data, raw_data[["gene_id", "reverse"]]], axis=1)
                    data.columns=['gene_id', basename]
                else:
                    raw_data.columns = ['gene_id', basename]
                    data = pd.merge(left = data, right = raw_data, how='left', on="gene_id")
    return data

data = combine_star_counts(star_dir)
data.to_csv('deseq2_input.csv')
