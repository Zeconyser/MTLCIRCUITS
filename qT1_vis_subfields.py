import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# --- Paths ---
base = "/Volumes/SFB/B04/MTLCIRCUITS/Tables"
csv_path = os.path.join(base,"MTLCIRCUITS_7T.csv")
out_dir = os.path.join(base, "layer_profile_plots")
os.makedirs(out_dir, exist_ok=True)

rois = ["ERC", "SUB", "CA1", "CA2", "CA3", "DG"]
hemis = ["left", "right"]
n_layers = 21

FLIP_LAYERS = False  # Set True only if inner/outer surfaces were swapped during rim generation

# --- Load wide-format CSV ---
df_wide = pd.read_csv(csv_path)
df_wide["ID"] = df_wide["ID"].astype(str)

# --- Assign Layer Zones (same thresholds used earlier) ---
def assign_zone(layer):
    if 1 <= layer <= 7:
        return "Inner"
    elif 8 <= layer <= 14:
        return "Middle"
    elif 15 <= layer <= 21:
        return "Outer"
    return None


def melt_roi_to_long(df_wide, roi):
    """
    Reshapes the wide-format layer columns for one ROI into long format:
    Subject, Hemisphere, Layer, Myelin
    """
    records = []
    for hemi in hemis:
        layer_cols = [f"mean_{roi}_layer{l}_{hemi}" for l in range(1, n_layers + 1)]
        present_cols = [c for c in layer_cols if c in df_wide.columns]

        for _, row in df_wide.iterrows():
            for col in present_cols:
                layer_num = int(col.split("layer")[1].split(f"_{hemi}")[0])
                value = row[col]
                if pd.notna(value):
                    records.append({
                        "Subject": row["ID"],
                        "Hemisphere": hemi,
                        "Layer": layer_num,
                        "Myelin": value
                    })
    return pd.DataFrame(records)


# --- Plot styling (matches original CA1 script) ---
layer_zones = {"Inner": (1, 7.5), "Middle": (7.5, 14.5), "Outer": (14.5, 21)}
colors = {"Inner": "#ffcccc", "Middle": "#ffff99", "Outer": "#add8e6"}

for roi in rois:
    df = melt_roi_to_long(df_wide, roi)

    if df.empty:
        print(f"No layer data found for {roi}, skipping.")
        continue

    if FLIP_LAYERS:
        df["Layer_flipped"] = (n_layers + 1) - df["Layer"]
    else:
        df["Layer_flipped"] = df["Layer"]

    df["LayerZone"] = df["Layer_flipped"].apply(assign_zone)
    df = df.dropna(subset=["LayerZone", "Myelin"])

    # --- Summary table ---
    summary = df.groupby("LayerZone")["Myelin"].agg(["mean", "std", "count"]).reset_index()
    summary = summary.rename(columns={"mean": "Mean_Myelin", "std": "STD_Myelin", "count": "N"})
    summary["LayerZone"] = pd.Categorical(summary["LayerZone"], categories=["Inner", "Middle", "Outer"], ordered=True)
    summary = summary.sort_values("LayerZone").rename(columns={"LayerZone": "Compartment"})
    print(f"\n--- {roi} summary ---")
    print(summary)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharey=True, sharex=True)

    for ax, hemi in zip(axes, hemis):
        hemi_df = df[df["Hemisphere"] == hemi]

        # Individual subjects
        for subject in hemi_df["Subject"].unique():
            subj_df = hemi_df[hemi_df["Subject"] == subject].sort_values("Layer_flipped")
            ax.plot(subj_df["Layer_flipped"], subj_df["Myelin"], color='black', alpha=0.5)

        # Zone shading
        for zone, (start, end) in layer_zones.items():
            ax.axvspan(start, end, color=colors[zone], alpha=0.5, zorder=0)

        # Group mean
        group_mean = hemi_df.groupby("Layer_flipped")["Myelin"].mean()
        ax.plot(group_mean.index, group_mean.values, color='red', lw=2, label='Group Mean')

        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=18)

        ax.set_title(f"qT1 Profile - {hemi.capitalize()} {roi}", fontweight='bold', size=24)
        ax.set_xlabel("Layer", size=18)
        ax.set_ylabel("qT1 (ms)", size=18)

        ax.set_xlim(1, n_layers)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    out_path = os.path.join(base, f"qt1_layer_profile_{roi}.png")
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_path}")

print("\nDone.")
