import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

mega_file = pd.DataFrame()

def find_quant_files(directory):
    quant_files = []
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith("_quant.sf"):
                quant_files.append(os.path.join(dirpath, filename))
    return quant_files

quant_files = find_quant_files("salmon_quantification")

for file in quant_files:
    sample_name = os.path.basename(file).split("_")[0]
    data = pd.read_csv(file, sep="\t", header=0)
    data["sample"] = sample_name
    mega_file = pd.concat([mega_file, data], axis=0)

top_file = mega_file.groupby("Name").sort_values("TPM", ascending=False)
top_file = top_file.head(1000)

g = sns.scatterplot(data = top_file, x = "Name", y = "TPM", hue = "sample")
g.set_yscale("log")
g = plt.savefig("top_1000_genes.png")
