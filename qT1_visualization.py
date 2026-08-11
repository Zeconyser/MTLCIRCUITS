import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
    First-look visualizations for subfield qT1 means:
      1. Boxplots of ROI mean qT1, one panel per hemisphere
      2. Correlation matrix between ROIs, one per hemisphere
'''

base = "/Volumes/SFB/B04/MTLCIRCUITS"
csv_path = os.path.join(base, "Tables",  "MTLCIRCUITS_7T.csv")

df = pd.read_csv(csv_path)

hemis = ["left", "right"]
rois = ["ERC", "SUB", "CA1", "CA2", "CA3", "DG"]

sns.set_theme(style="whitegrid", context="talk")

# ---------------------------------------------------------
# 1. Boxplots of ROI means, one subplot per hemisphere
# ---------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, hemi in zip(axes, hemis):
    roi_cols = [f"mean_{roi}_{hemi}" for roi in rois]
    present_cols = [c for c in roi_cols if c in df.columns]

    plot_df = df[present_cols].melt(var_name="ROI", value_name="qT1")
    plot_df["ROI"] = plot_df["ROI"].str.replace(f"mean_", "").str.replace(f"_{hemi}", "")

    sns.boxplot(data=plot_df, x="ROI", y="qT1", ax=ax, order=rois)
    sns.stripplot(data=plot_df, x="ROI", y="qT1", ax=ax, order=rois,
                  color="black", alpha=0.4, size=3, jitter=True)

    ax.set_title(f"{hemi.capitalize()} hemisphere")
    ax.set_xlabel("")
    ax.set_ylabel("Mean qT1" if hemi == "left" else "")

fig.suptitle("Subfield mean qT1 by hemisphere", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(base, "roi_boxplots.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------
# 2. Correlation matrices, one per hemisphere
# ---------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, hemi in zip(axes, hemis):
    roi_cols = [f"mean_{roi}_{hemi}" for roi in rois]
    present_cols = [c for c in roi_cols if c in df.columns]

    corr = df[present_cols].corr()
    corr.columns = [c.replace("mean_", "").replace(f"_{hemi}", "") for c in corr.columns]
    corr.index = corr.columns

    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", vmin=-1, vmax=1,
                square=True, ax=ax, cbar_kws={"shrink": 0.8})
    ax.set_title(f"{hemi.capitalize()} hemisphere")

fig.suptitle("Subfield qT1 correlation matrix", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(base,"Tables", "roi_correlation_matrices.png"), dpi=200, bbox_inches="tight")
plt.close(fig)

print("Saved:")
print(f"  {os.path.join(base, 'roi_boxplots.png')}")
print(f"  {os.path.join(base, 'roi_correlation_matrices.png')}")
