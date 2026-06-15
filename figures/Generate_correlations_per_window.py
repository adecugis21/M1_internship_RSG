#!/usr/bin/env python3
import os
import sys
import re
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import pyBigWig
from tqdm import tqdm

# ------------------------------------------------------------------
# Import Momics
# ------------------------------------------------------------------
sys.path.insert(0, "/pasteur/appa/scratch/adecugis/all_momics/momics/src")
from momics.momics import Momics


"""

python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/correlations/fun_graphs_MNASE.py/V3_make_MNASE_Manon_big_graph_corr_CONCATENATED_DATA.py \
  --pred-path /pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/fun_graph/Pombe_modulable_June10_Manon_3_species_no_maxpool \
  --summary-path /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/test_region_summary.txt \
  --n-windows-species 500


"""


# ------------------------------------------------------------------
# Species / repo mapping
# ------------------------------------------------------------------
SPECIES_REPO = {
    "S_pombe":      "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/S_pombe.momics",
    "S_cerevisiae": "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/S_cerevisiae.momics",
    "C_glabrata":   "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/c_glabrata.momics",
    "S_eubayanus":  "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/S_eubayanus.momics",
    "S_paradoxus":  "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/S_paradoxus.momics",
}

PREFIX_MAP = {
    "S_pombe": ("972hminus", "972hplus", "X54421"),
    "C_glabrata": ("CBS138", "C_glabrata"),
    "S_cerevisiae": ("S288c", "S_cerevisiae"),
    "S_eubayanus": ("CDS216", "S_eubayanus"),
    "S_paradoxus": ("CBS432", "S_paradoxus",),
}

# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument(
    "--pred-path",
    type=str,
    required=True,
    help="Root folder containing prediction BigWigs (searched recursively)",
)
parser.add_argument(
    "--summary-path",
    type=str,
    default=None,
    help="Path to test_region_summary.txt (auto-discovered if omitted)",
)
parser.add_argument(
    "--n-windows-species",
    type=int,
    default=5000,
    help="Number of windows to sample evenly across test regions for each species",
)
args = parser.parse_args()

PRED_ROOT = Path(args.pred_path).resolve()

def find_summary_path(pred_root: Path, explicit: str | None) -> Path:
    if explicit:
        return Path(explicit).resolve()
    for p in [pred_root, *pred_root.parents]:
        candidate = p / "test_region_summary.txt"
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not find test_region_summary.txt; pass --summary-path")

SUMMARY_PATH = find_summary_path(PRED_ROOT, args.summary_path)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

import warnings
from scipy.stats import ConstantInputWarning  # available scipy ≥ 1.9

# ------------------------------------------------------------------
# Helper: safe correlation computation
# ------------------------------------------------------------------
def safe_correlations(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    """
    Compute Pearson and Spearman r.
    Returns (nan, nan) if either array is constant, with a clear reason logged.
    """
    if np.std(y_true) == 0 or np.std(y_pred) == 0:
        # Identify which side is constant for better diagnostics
        which = []
        if np.std(y_true) == 0:
            which.append(f"y_true (constant={y_true[0]:.4g})")
        if np.std(y_pred) == 0:
            which.append(f"y_pred (constant={y_pred[0]:.4g})")
        print(f"    ⚠️  Constant array detected: {', '.join(which)} — skipping correlation")
        return np.nan, np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("error", ConstantInputWarning)
        try:
            pear_r, _ = pearsonr(y_true, y_pred)
            spear_r, _ = spearmanr(y_true, y_pred)
        except ConstantInputWarning as exc:
            print(f"    ⚠️  ConstantInputWarning caught: {exc}")
            return np.nan, np.nan

    return pear_r, spear_r

def infer_species_from_chrom(chrom: str) -> str | None:
    for species, prefixes in PREFIX_MAP.items():
        if chrom.startswith(prefixes):
            return species
    return None


def parse_test_regions(summary_path: Path) -> list[dict]:
    text = summary_path.read_text()
    regions = []
    repo_blocks = re.split(r"=== Test \(", text)[1:]

    for block in repo_blocks:
        repo_match = re.match(r"(.+?)\) bin ranges by chromosome ===", block)
        if not repo_match:
            continue

        for line in block.split("\n"):
            if "bins_count" in line or "===" in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                chrom = parts[0]
                try:
                    regions.append(
                        {
                            "chrom": chrom,
                            "start": int(parts[2]),
                            "end": int(parts[3]),
                            "species": infer_species_from_chrom(chrom),
                        }
                    )
                except ValueError:
                    continue
    return regions


def group_regions_by_species(regions: list[dict]) -> dict[str, list[dict]]:
    grouped = {}
    for r in regions:
        species = r["species"]
        if species is None:
            continue
        grouped.setdefault(species, []).append(r)
    return grouped


def build_species_window_plan(
    regions: list[dict],
    n_windows: int,
    window_size: int = 512,
    stride: int = 48,
) -> list[tuple[str, int, int]]:
    windows = []
    for r in regions:
        start, end = r["start"], r["end"]
        if (end - start) < window_size:
            continue
        starts = np.arange(start, end - window_size + 1, stride, dtype=int)
        windows.extend((r["chrom"], int(s), int(s + window_size)) for s in starts)

    total = len(windows)
    if total == 0:
        raise ValueError("No valid windows could be built from test regions.")
    if n_windows > total:
        raise ValueError(
            f"Requested --n-windows-species={n_windows}, but only {total} windows available."
        )

    # evenly spread indices over the full window list (unique, deterministic)
    idx = (np.arange(n_windows) * total) // n_windows
    return [windows[i] for i in idx]


def discover_prediction_files(pred_root: Path) -> dict[str, dict[str, Path]]:
    """
    Expected filenames:
      pred_<anything>_<target>_softplus.bw
    Works with patterns like:
      pred_CBS432_chrXVI_MNASE_15_softplus.bw  -> MNASE_15
      pred_972hminus_ENA|CU329672|CU329672.1_H3_softplus.bw  -> H3
    """
    grouped = {}

    for bw_file in pred_root.rglob("*.bw"):
        # Infer species name from folder path (case-insensitive)
        species = next((sp for sp in SPECIES_REPO if sp.lower() in [p.lower() for p in bw_file.parts]), None)
        if species and species not in grouped:
            grouped[SPECIES_REPO.get(species, species)] = {}
        if species is None:
            continue

        # -----------------------------------------------------------
        # 1️⃣ match any prefix pred_*_XXX_softplus.bw
        # -----------------------------------------------------------
        m = re.match(r"^pred_.+_([A-Za-z0-9_]+)_softplus\.bw$", bw_file.name)
        if not m:
            m = re.match(r"^pred_.+_([A-Za-z0-9_]+)\.bw$", bw_file.name)
        if not m:
            continue

        raw_target = m.group(1)

        # -----------------------------------------------------------
        # 2️⃣ Normalize: fix MNASE_<num> pattern
        # -----------------------------------------------------------
        # If it looks like "<number>", add MNASE_ prefix
        if re.fullmatch(r"\d+", raw_target):  # 15, 25, 35
            target = f"MNASE_{raw_target}"
        # If it ends with "_<number>" but missing MNASE, prepend MNASE_
        elif re.fullmatch(r"[A-Za-z]+_(\d+)", raw_target) and not raw_target.startswith("MNASE"):
            target = f"MNASE_{raw_target.split('_')[-1]}"
        # Keep if it already has MNASE_X
        elif raw_target.startswith("MNASE_"):
            target = raw_target
        else:
            target = raw_target  # e.g. H3, H3K4ME3

        print(f"[{species}] Found target: {target}")
        grouped.setdefault(species, {})[target] = bw_file

    return grouped


def load_bigwig_dict(bw_file: Path) -> dict[str, np.ndarray]:
    with pyBigWig.open(str(bw_file)) as bw:
        return {
            chrom: np.array(bw.values(chrom, 0, length))
            for chrom, length in bw.chroms().items()
        }


def get_signal(repo, track, chrom, start, end):
    covs = repo.tracks(track)
    arr = np.array(covs[chrom])
    seg = arr[start:end]
    return seg[np.isfinite(seg)]


# ------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------
test_regions = parse_test_regions(SUMMARY_PATH)
species_regions = group_regions_by_species(test_regions)
pred_files = discover_prediction_files(PRED_ROOT)

print("SUMMARY ----------------------------")
print(f"Prediction root : {PRED_ROOT}")
print(f"Summary path     : {SUMMARY_PATH}")
print(f"Species found    : {sorted(species_regions.keys())}")

species = sorted(species_regions.keys())
# ------------------------------------------------------------------
# Compute correlations
# ------------------------------------------------------------------
species_target_records = []


for species, regions in species_regions.items():
    if species not in SPECIES_REPO:
        continue
    if species not in pred_files:
        print(f"⚠️ No prediction files found for {species}")
        continue

    repo = Momics(SPECIES_REPO[species])
    window_plan = build_species_window_plan(regions, args.n_windows_species)
    print(
        f"\nSpecies: {species} | regions={len(regions)} | "
        f"sampled_windows={len(window_plan)} | targets={len(pred_files[species])}"
    )
    print(window_plan)

    for target, bw_file in sorted(pred_files[species].items()):
        print(f"  Loading {target} -> {bw_file.name}")
        print(f"  full path -> {bw_file}")
        predictions = load_bigwig_dict(bw_file)

        y_pred_all, y_true_all = [], []
        pbar = tqdm(total=len(window_plan), desc=f"{species}:{target}", leave=False)

        for chrom, s, e in window_plan:
            if chrom not in predictions:
                pbar.update(1)
                continue

            pred_chrom = np.array(predictions[chrom])
            if e > len(pred_chrom):
                pbar.update(1)
                continue

            y_pred_all.append(np.nanmax(pred_chrom[s:e]))
            arr_true = get_signal(repo, target, chrom, s, e)
            y_true_all.append(np.nanmax(arr_true) if arr_true.size else np.nan)
            pbar.update(1)

        pbar.close()

        y_pred_all = np.array(y_pred_all)
        y_true_all = np.array(y_true_all)
        mask = np.isfinite(y_true_all) & np.isfinite(y_pred_all)
        y_pred_all = y_pred_all[mask]
        y_true_all = y_true_all[mask]

        if len(y_true_all) < 2:
            print(f"⚠️ Skipping {species}:{target} because no valid overlapping windows")
            continue

        pear_r, spear_r = safe_correlations(y_true_all, y_pred_all)

        if np.isnan(pear_r) and np.isnan(spear_r):
            print(f"  ⚠️  Skipping record for {species}:{target} (all-NaN correlations)")
            continue  # or keep the record with NaN — your choice

        species_target_records.append(
            {
                "species": species,
                "target": target,
                "pearson_r": pear_r,
                "spearman_r": spear_r,
                "n_windows": len(y_true_all),
            }
        )

# ------------------------------------------------------------------
# Save summary tables
# ------------------------------------------------------------------
run_name = PRED_ROOT.name
run_with_step = f"{run_name}_{args.n_windows_species}_windows"
window_number = f"{args.n_windows_species}_windows"

out_dir = Path("/pasteur/helix/scratch/adecugis/all_momics/momics/output/figures") / "fun_figure" / window_number / run_with_step
out_dir.mkdir(parents=True, exist_ok=True)

by_target = pd.DataFrame(species_target_records)
by_target.to_csv(out_dir / "summary_correlation_by_target.csv", index=False)

species_summary = (
    by_target.groupby("species", as_index=False)
    .agg(
        pearson_r=("pearson_r", "mean"),
        spearman_r=("spearman_r", "mean"),
        n_windows=("n_windows", "sum"),
        n_targets=("target", "count"),
    )
)

species_summary.to_csv(out_dir / "summary_correlation_by_species.csv", index=False)

print(f"✔ Saved: {out_dir / 'summary_correlation_by_target.csv'}")
print(f"✔ Saved: {out_dir / 'summary_correlation_by_species.csv'}")
