import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# no_cleanup = pd.read_csv("multiqc_test/Cappable-seq-pilot_multiqc_report_data/multiqc_salmon_salmon_no_cleanup.txt", 
#                          sep="\t")
# with_cleanup = pd.read_csv("multiqc_test/Cappable-seq-pilot_multiqc_report_data/multiqc_salmon_salmon_after_cleanup.txt",
#                          sep="\t")

# trinity = pd.read_csv("multiqc_test/Cappable-seq-pilot_multiqc_report_data/multiqc_salmon_salmon_mapped_to_trinity.txt",
#                          sep="\t")

star = pd.read_csv("multiqc_test/Cappable-seq-pilot_multiqc_report_data/multiqc_star_star.txt",
                            sep="\t")

salmon_gencode = pd.read_csv("multiqc_test/Cappable-seq-pilot_multiqc_report_data/multiqc_salmon_salmon_gencode.txt",
                                sep="\t")


def graph(data, title):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    data['percent_mapped'] = data['num_mapped'] / data['num_processed'] * 100
    ax = sns.barplot(x="Sample", y="percent_mapped", data=data)
    ax.set_title(f"{title} Reads Mapped by Salmon")
    ax.set_ylabel("Percent of Reads Mapped")
    ax.set_xlabel("Sample")
    plt.xticks(rotation=45)
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{title}_mapped.png", dpi=300)


# graph(no_cleanup, "No Cleanup")
# graph(with_cleanup, "With Cleanup")
# graph(trinity, "Trinity")
graph(salmon_gencode, "Salmon Gencode")

# Include STAR as one off
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 6))
ax = sns.barplot(x="Sample", y="uniquely_mapped_percent", data=star)
ax.set_title("STAR Reads Aligned to genome")
ax.set_ylabel("Percent of Reads Aligned")
ax.set_xlabel("Sample")
plt.xticks(rotation=45)
plt.ylim(0, 100)
plt.tight_layout()
plt.savefig("STAR_mapped.png", dpi=300)