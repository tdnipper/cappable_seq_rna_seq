import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load the data
data_infection = pd.read_csv("geneCounts_infection.csv")
filtered_infection = data_infection[data_infection["padj"] <= 0.05]

# Plot the data
sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Plot infection

ax1 = axes[0]
sns.scatterplot(
    data=data_infection,
    x="log2FoldChange",
    y="baseMean",
    size="lfcSE",
    linewidth=0,
    legend=False,
    ax=ax1
)
ax1.set_title("Raw Data")

ax2 = axes[1]
sns.scatterplot(
    data=filtered_infection,
    x="log2FoldChange",
    y="baseMean",
    size="lfcSE",
    linewidth=0,
    legend=False,
    ax=ax2
)
ax2.set_title("Filtered data")
plt.xlabel("log2FoldChange (p<=0.05)")
fig.suptitle("Fold change vs abundance during infection and enrichment")
plt.tight_layout()
plt.savefig("infection_test.png", dpi=300)
