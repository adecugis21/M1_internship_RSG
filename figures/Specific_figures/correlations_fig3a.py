#!/usr/bin/env python3
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import pyBigWig

# ------------------------------------------------------------------
# Fixed setup: use your BigWigs directly
# ------------------------------------------------------------------
species = "cerevisiae"
TARGET = "H3K4ME3"
STEP_SAMPLE = 50
TARGET_SIZE = 512
STRIDE = 48

bw_folder = Path(
    "/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/2_species_3tracks_june2"
) / species
bw_folder_single = Path(
    "/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/2_species_1_track_june2"
) / species
gt_bw = Path(
    "/pasteur/helix/scratch/adecugis/all_momics/momics/output/a_extracting_bw"
) / species / f"{species}_{TARGET}_rescaled.bw"
pred_single_bw = bw_folder_single / f"{species}_pred_XV_{TARGET}.bw"
pred_multi_bw = bw_folder / f"{species}_pred_XV_{TARGET}.bw"

OUTPUT_DIR = (
    Path("/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/correlation/figures")
    / f"{species}_corr_bw_only"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
corr_dir = OUTPUT_DIR / "window_corr"
corr_dir.mkdir(exist_ok=True)
graph_dir = corr_dir / "graphs"
graph_dir.mkdir(exist_ok=True)

print(f"GT   : {gt_bw}")
print(f"single: {pred_single_bw}")
print(f"multi : {pred_multi_bw}")
print(f"Output: {OUTPUT_DIR}")

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def load_bw(bw_path: Path):
    with pyBigWig.open(str(bw_path)) as bw:
        return {chrom: np.array(bw.values(chrom, 0, length)) for chrom, length in bw.chroms().items()}

def window_correlations(gt_chroms, pred_chroms):
    records = []
    plot_data = []

    shared_chroms = sorted(set(gt_chroms) & set(pred_chroms))
    shared_chroms = [c for c in shared_chroms if "XV" in c]  # Only XV
    
    for chrom in shared_chroms:
        gt = gt_chroms[chrom]
        pred = pred_chroms[chrom]
        max_len = min(len(gt), len(pred))
        starts = np.arange(0, max_len - TARGET_SIZE + 1, STRIDE)
        starts = starts[::STEP_SAMPLE]

        y_true, y_pred = [], []
        for s in starts:
            e = s + TARGET_SIZE
            if e > max_len:
                break
            gt_win = gt[s:e]
            pred_win = pred[s:e]

            gt_win = gt_win[np.isfinite(gt_win)]
            pred_win = pred_win[np.isfinite(pred_win)]
            if gt_win.size == 0 or pred_win.size == 0:
                continue

            y_true.append(np.nanmax(gt_win))
            y_pred.append(np.nanmax(pred_win))

        if len(y_true) < 2:
            continue

        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        mask = np.isfinite(y_true) & np.isfinite(y_pred)
        y_true = y_true[mask]
        y_pred = y_pred[mask]
        if len(y_true) < 2:
            continue

        pear_r, _ = pearsonr(y_true, y_pred)
        spear_r, _ = spearmanr(y_true, y_pred)

        records.append(
            {
                "chrom": chrom,
                "pearson_r": pear_r,
                "spearman_r": spear_r,
                "n_windows": len(y_true),
            }
        )
        plot_data.append((chrom, y_true, y_pred, pear_r, spear_r))

        pd.DataFrame({"y_true_max": y_true, "y_pred_max": y_pred}).to_csv(
            corr_dir / f"{chrom}_max_per_window.csv", index=False
        )

    return records, plot_data

def save_scatter_plot(y_true, y_pred, title, outpath):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(y_true, y_pred, s=10, alpha=0.6)
    # x-axis stays 0..1; set y-axis from preds so small values are visible
    xlim = 1.0
    # use 99.9th percentile of predictions, with small minimum and cap at 1.0
    try:
        p99 = float(np.nanpercentile(y_pred, 99.9))
    except Exception:
        p99 = 0.0
    ymax = max(1e-3, p99 * 1.2)  # ensure tiny predictions are visible
    ymax = min(1.0, ymax)        # keep within [0,1]
    ax.plot([0, xlim], [0, xlim], color="red", linewidth=1)
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("True max (coverage per 512 bp)")
    ax.set_ylabel("Pred max (coverage per 512 bp)")
    ax.set_title(title)
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close(fig)

# ------------------------------------------------------------------
# Load BigWigs
# ------------------------------------------------------------------
gt = load_bw(gt_bw)
pred_single = load_bw(pred_single_bw)
pred_multi = load_bw(pred_multi_bw)

# ------------------------------------------------------------------
# Correlations
# ------------------------------------------------------------------
summary = []
plot_data_by_label = {}

for label, pred in [
    ("single_track", pred_single),
    ("multi_modal", pred_multi),
]:
    records, plot_data = window_correlations(gt, pred)
    # convert list -> dict for later combined plotting
    plot_dict = {chrom: (y_true, y_pred, pear_r, spear_r) for chrom, y_true, y_pred, pear_r, spear_r in plot_data}
    plot_data_by_label[label] = plot_dict

    if records:
        df = pd.DataFrame(records)
        df.insert(0, "comparison", label)
        summary.append(df)

        for chrom, y_true, y_pred, pear_r, spear_r in plot_data:
            save_scatter_plot(
                y_true,
                y_pred,
                f"{TARGET} - {label} - {chrom}\nPearson r={pear_r:.3f} | Spearman r={spear_r:.3f}",
                graph_dir / f"{TARGET}_{label}_{chrom}_scatter.png",
            )

        pd.DataFrame(records).to_csv(corr_dir / f"{label}_summary_correlation.csv", index=False)
        print(f"✔ Saved {label} outputs in {corr_dir}")
    else:
        print(f"⚠️ No valid windows for {label}")

# Combined plots per chromosome where both labels have data
labels = list(plot_data_by_label.keys())
if len(labels) >= 2:
    lab1, lab2 = labels[0], labels[1]
    colors = {lab1: "C0", lab2: "C1"}
    markers = {lab1: "o", lab2: "s"}
    shared_chroms = sorted(set(plot_data_by_label[lab1].keys()) & set(plot_data_by_label[lab2].keys()))
    for chrom in shared_chroms:
        y1, p1, pear1, spear1 = plot_data_by_label[lab1][chrom][:4]
        y2, p2, pear2, spear2 = plot_data_by_label[lab2][chrom][:4]
        y1 = np.asarray(y1); p1 = np.asarray(p1)
        y2 = np.asarray(y2); p2 = np.asarray(p2)
 
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(y1, p1, s=10, alpha=0.6, color=colors[lab1], marker=markers[lab1], label=f"{lab1} (r={pear1:.3f})")
        ax.scatter(y2, p2, s=10, alpha=0.6, color=colors[lab2], marker=markers[lab2], label=f"{lab2} (r={pear2:.3f})")
 
        xlim = 1.0
        # compute combined ymax from both prediction sets (99.9th percentile)
        try:
            p99_1 = float(np.nanpercentile(p1, 99.9))
        except Exception:
            p99_1 = 0.0
        try:
            p99_2 = float(np.nanpercentile(p2, 99.9))
        except Exception:
            p99_2 = 0.0
        ymax = max(1e-3, max(p99_1, p99_2) * 1.2)
        ymax = min(1.0, ymax)
 
         # draw best-fit line per model
        for lbl, yy, pp in [(lab1, y1, p1), (lab2, y2, p2)]:
             mask = np.isfinite(yy) & np.isfinite(pp)
             if mask.sum() >= 2:
                slope, intercept = np.polyfit(yy[mask], pp[mask], 1)
                xs = np.array([0, xlim])
                ax.plot(xs, slope * xs + intercept, color=colors[lbl], linewidth=1.5)
 
        ax.plot([0, xlim], [0, xlim], color="gray", linewidth=0.8, linestyle="--")
        ax.set_xlim(0, xlim)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("True max (coverage per 512 bp)")
        ax.set_ylabel("Pred max (coverage per 512 bp)")
        ax.set_title(f"{TARGET} - combined - {chrom}")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.savefig(graph_dir / f"{TARGET}_combined_{chrom}_scatter.png", dpi=300)
        plt.close(fig)

if summary:
    pd.concat(summary, ignore_index=True).to_csv(corr_dir / "summary_correlation.csv", index=False)

print(f"✔ Finished. Plots and CSVs saved to {graph_dir}")