#!/usr/bin/env python3
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import matplotlib as mpl


# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
BASE_DIR = Path("/pasteur/helix/scratch/adecugis/all_momics/momics/output/figures/fun_figure/200_windows")

# ------------------------------------------------------------------
# Helper: infer species name from run_name
# ------------------------------------------------------------------
def infer_species(run_name: str) -> str:
    """
    Guess species name based on common patterns in the run_name.
    You can expand this list depending on your naming conventions.
    """
    run_name_lower = run_name.lower()
    if "cerevis" in run_name_lower or "s288" in run_name_lower:
        return "S_cerevisiae"
    elif "pombe" in run_name_lower or "972h" in run_name_lower:
        return "S_pombe"
    elif "glabrata" in run_name_lower or "cbs138" in run_name_lower:
        return "C_glabrata"
    elif "paradoxus" in run_name_lower:
        return "S_paradoxus"
    elif "eubayanus" in run_name_lower:
        return "S_eubayanus"
    elif "concatenated" in run_name_lower:
        return "concatenated"
    else:
        return "unknown"

# ------------------------------------------------------------------
# Global font sizes
# ------------------------------------------------------------------
FONT_SIZE        = 19   # base size — change this one number to scale everything
TITLE_SIZE       = FONT_SIZE + 4
LABEL_SIZE       = FONT_SIZE + 2
TICK_SIZE        = FONT_SIZE
ANNOT_SIZE       = FONT_SIZE - 2   # numbers inside heatmap cells
CBAR_LABEL_SIZE  = FONT_SIZE

plt.rcParams.update({
    "font.size":        FONT_SIZE,
    "axes.titlesize":   TITLE_SIZE,
    "axes.labelsize":   LABEL_SIZE,
    "xtick.labelsize":  TICK_SIZE,
    "ytick.labelsize":  TICK_SIZE,
})

# ------------------------------------------------------------------
# Load all correlation data from subfolders
# ------------------------------------------------------------------
all_records = []

for subfolder in sorted(BASE_DIR.iterdir()):
    if not subfolder.is_dir():
        continue
    
    csv_file = subfolder / "summary_correlation_by_target.csv"
    if not csv_file.exists():
        print(f"⚠️ Skipping {subfolder.name} (no summary_correlation_by_target.csv)")
        continue
    
    try:
        df = pd.read_csv(csv_file)
        df["run_name"] = subfolder.name
        df["species_label"] = df["run_name"].apply(infer_species)  # 🧬 apply function here
        all_records.append(df)
        print(f"✔ Loaded: {subfolder.name}")
    except Exception as e:
        print(f"✗ Error loading {subfolder.name}: {e}")

if not all_records:
    print("No correlation files found!")
    sys.exit(1)

combined_df = pd.concat(all_records, ignore_index=True)
print(f"\nTotal records: {len(combined_df)}")
print(f"Unique species in CSV: {combined_df['species'].unique()}")
print(f"Inferred from run names: {combined_df['species_label'].unique()}")
print(f"Unique targets: {combined_df['target'].unique()}")

# Create output directory
output_dir = BASE_DIR / "output_graphs"
output_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Create per-target heatmaps: runs x species
# ------------------------------------------------------------------
for target in sorted(combined_df["target"].unique()):
    target_df = combined_df[combined_df["target"] == target].copy()

    # make species_label an ordered categorical with "concatenated" first
    order = ["concatenated", "S_pombe", "S_paradoxus", "S_eubayanus",
             "S_cerevisiae", "C_glabrata"]
    target_df["species_label"] = pd.Categorical(
        target_df["species_label"], categories=order, ordered=True
    )
    
    fig, axes = plt.subplots(1, 2, figsize=(18, 9))   # wider to fit bigger text
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)   # ← add this line

    # Pearson heatmap
    pearson_pivot = target_df.pivot_table(
        values="pearson_r",
        index="species_label",
        columns="species",
        aggfunc="mean"
    ).sort_index(level=0)

    sns.heatmap(
        pearson_pivot,
        annot      = True,
        fmt        = ".2f",
        annot_kws  = {"size": ANNOT_SIZE},   # ← cell numbers
        cmap       = "RdYlGn",
        center     = 0,
        ax         = axes[0],
        vmin       = -1,
        vmax       = 1,
        cbar_kws   = {"label": "Pearson r", "shrink": 0.8}
    )
    axes[0].set_title(f"{target} — Pearson Correlation", fontsize=TITLE_SIZE, pad=12)
    axes[0].set_xlabel("Species of the test data",       fontsize=LABEL_SIZE)
    axes[0].set_ylabel("Model training data",            fontsize=LABEL_SIZE)
    axes[0].tick_params(axis="both", labelsize=TICK_SIZE)

    # Fix colorbar label size
    axes[0].collections[0].colorbar.ax.tick_params(labelsize=TICK_SIZE)
    axes[0].collections[0].colorbar.set_label("Pearson r", fontsize=CBAR_LABEL_SIZE)

    # Spearman heatmap
    spearman_pivot = target_df.pivot_table(
        values="spearman_r",
        index="species_label",
        columns="species",
        aggfunc="mean"
    ).sort_index(level=0)

    sns.heatmap(
        spearman_pivot,
        annot      = True,
        fmt        = ".2f",
        annot_kws  = {"size": ANNOT_SIZE},
        cmap       = "RdYlGn",
        center     = 0,
        ax         = axes[1],
        vmin       = -1,
        vmax       = 1,
        cbar_kws   = {"label": "Spearman r", "shrink": 0.8}
    )
    axes[1].set_title(f"{target} — Spearman Correlation", fontsize=TITLE_SIZE, pad=12)
    axes[1].set_xlabel("Species of the test data",        fontsize=LABEL_SIZE)
    axes[1].set_ylabel("Model training data",             fontsize=LABEL_SIZE)
    axes[1].tick_params(axis="both", labelsize=TICK_SIZE)

    """
            # Add horizontal colorbar below axes[1]
    sm2 = mpl.cm.ScalarMappable(cmap="RdYlGn", norm=norm)
    sm2.set_array([])
    cbar1 = fig.colorbar(sm2, ax=axes[1], orientation="horizontal",
                        fraction=0.046, pad=2, aspect=30)
    cbar1.set_label("Spearman r", fontsize=CBAR_LABEL_SIZE)
    cbar1.ax.tick_params(labelsize=TICK_SIZE)
    """

    axes[1].collections[0].colorbar.ax.tick_params(labelsize=TICK_SIZE)
    axes[1].collections[0].colorbar.set_label("Spearman r", fontsize=CBAR_LABEL_SIZE)

    plt.tight_layout()
    out_path = output_dir / f"heatmap_{target}_by_inferred_species.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"✔ Saved: {out_path}")
    plt.close()

print("\n✔ All visualizations complete!")