import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load the data
data_condition = pd.read_csv("res05_condition_annotated.csv")
data_infection = pd.read_csv("res05_infection_annotated.csv")

# Plot the data
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Plot for data_condition
ax1 = axes[0]
sns.scatterplot(
    data=data_condition,
    x="log2FoldChange",
    y="baseMean",
    size="lfcSE",
    ax=ax1,
    linewidth=0,
    legend=False,
)
ax1.set_title("Cappable-seq recovery")

# Plot for data_infection
ax2 = axes[1]
sns.scatterplot(
    data=data_infection,
    x="log2FoldChange",
    y="baseMean",
    size="lfcSE",
    ax=ax2,
    linewidth=0,
    legend=False,
)
ax2.set_title("Infection")

plt.tight_layout()
plt.savefig("fc_vs_abundance.png", dpi=300)
plt.close()
# Plot data for both log2FC
fig2 = sns.regplot(
    x=data_condition["log2FoldChange"],
    y=data_infection["log2FoldChange"],
)
fig2.set(xlabel="Cappable-seq log2FoldChange", ylabel="Infection log2FoldChange")
fig2.set_title("Log2FoldChange comparison Capable-seq vs Infection")
plt.savefig("log2fc_cappable_seq_vs_infection.png", dpi=300)
