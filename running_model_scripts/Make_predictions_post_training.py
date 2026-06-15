#!/usr/bin/env python3
"""
Prediction script driven by a training run folder name.

Example:
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/modulable_concatenate/using_modulable_concatenated_genome.py \
  --run_name Pombe_modulable_2MNASE_Manon_and_lab_june4_no_maxpool \
  --base_train_dir /pasteur/helix/scratch/adecugis/all_momics/momics/output \
  --base_pred_dir  /pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions/figure_3b

  """

# ============================================================
# Imports
# ============================================================
import os
import sys
import re
import gc
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import pyranges as pr

import tensorflow as tf

sys.path.insert(0, "/pasteur/appa/scratch/adecugis/all_momics/momics/src")

from momics.momics import Momics
from momics.streamer import MomicsStreamer
from momics import aggregate as mma
import momics.nn_changed_window as nn_changed

# ============================================================
# Helper: match target name to prediction key (handles suffixes like _softplus)
# ============================================================
def find_prediction_key(target: str, predictions: dict) -> str | None:
    """
    Find a key in predictions dict that matches target name.
    Handles suffixes in either direction:
    - target='BDF2_softplus', predictions has 'BDF2' → match
    - target='BDF2', predictions has 'BDF2_softplus' → match
    """
    # 1. exact match
    if target in predictions:
        return target

    # 2. target has suffix, predictions has clean name
    # strip known suffixes from target and try again
    clean_target = re.sub(r'_(softplus|head|output)$', '', target)
    if clean_target in predictions:
        return clean_target

    # 3. predictions have suffix, target is clean
    matches = [k for k in predictions if target in k]
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        print(f"  ⚠️ Multiple keys match '{target}': {matches} — using first")
        return matches[0]

    return None

def detect_species(chrom_name: str) -> str:
    """
    Map a chromosome name to a species ID.
    Tries substring match first, then falls back to exact match.
    Returns 'unknown' if no match found.
    """
    # 1. Try substring patterns first (longer, more specific)
    for pattern, species in SPECIES_PATTERNS_SUBSTR.items():
        if pattern in chrom_name:
            return species

    # 2. Try exact match for short/generic chromosome names
    if chrom_name in SPECIES_PATTERNS_EXACT:
        return SPECIES_PATTERNS_EXACT[chrom_name]

    return "unknown"

# ============================================================
# Species patterns
# chromosome name substrings → species/repo lookup
# ============================================================
SPECIES_PATTERNS_SUBSTR = {
    "CBS138":  "C_glabrata",
    "S288C":   "S_cerevisiae",
    "S288c":   "S_cerevisiae",
    "CDS216":  "S_eubayanus",
    "CBS432":  "S_paradoxus",
    "972h":    "S_pombe",
    "CU329":   "S_pombe",     # ENA|CU329670|… names
}

# Exact match patterns (too short/generic for substring matching)
SPECIES_PATTERNS_EXACT = {
    "IV":  "S_cerevisiae",
    "A":   "C_glabrata",
}

# Map species id → momics repository path
SPECIES_REPO = {
    "S_pombe": "/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/pombe.momics",
    "S_cerevisiae": "/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/cerevisiae.momics",
    "C_glabrata": "/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/glabrata.momics",
    "concatenated": "/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/concatenated.momics",
    # add others as needed
}

# ============================================================
# Argument parsing
# ============================================================
parser = argparse.ArgumentParser(
    description="Run predictions on the test region from a training run.")
parser.add_argument("--run_name", type=str, required=True,
                    help="Name of training run subfolder, e.g. 'Manon_data_default_no_maxpool'")
parser.add_argument("--base_train_dir", type=str,
                    default="/pasteur/helix/scratch/adecugis/all_momics/momics/output",
                    help="Base directory where training run folders live")
parser.add_argument("--base_pred_dir", type=str,
                    default="/pasteur/helix/scratch/adecugis/all_momics/momics/output/predictions",
                    help="Base directory where prediction outputs will be saved")
parser.add_argument("--checkpoint", type=str, default=None,
                    help="Optional: explicit .keras checkpoint path. "
                         "If not given, the latest epoch checkpoint is used.")
parser.add_argument("--stride", type=int, default=48)
parser.add_argument("--batch_size", type=int, default=4)

args = parser.parse_args()

# ============================================================
# Locate run directory + checkpoint
# ============================================================
RUN_DIR = Path(args.base_train_dir) / args.run_name
CKPT_DIR = RUN_DIR / "checkpoints_stage1"
SUMMARY_FILE = RUN_DIR / "test_region_summary.txt"

if not RUN_DIR.exists():
    raise FileNotFoundError(f"Run directory not found: {RUN_DIR}")

if not SUMMARY_FILE.exists():
    raise FileNotFoundError(f"Test region summary not found: {SUMMARY_FILE}")

# Find checkpoint automatically if not supplied
if args.checkpoint:
    MODEL_PATH = Path(args.checkpoint)
else:
    all_ckpts = sorted(CKPT_DIR.glob("epoch_*.keras"))
    if not all_ckpts:
        raise FileNotFoundError(f"No checkpoints found in {CKPT_DIR}")
    MODEL_PATH = all_ckpts[-1]   # latest epoch
    print(f"Auto‑selected checkpoint: {MODEL_PATH}")

# ============================================================
# Parse test_region_summary.txt
# ============================================================
# Expected format (written by training script):
#
#  === Test regions (remaining bins) ===
#  /path/to/species.momics: 1234 test bins
#  === Test (/path/to/species.momics) bin ranges by chromosome ===
#  Chromosome                bins_count  min_start  max_end
#  ENA|CU329672|CU329672.1       1234     1200000  3456789
#

def parse_test_regions(summary_path):
    """
    Returns a list of dicts:
    [{"repo": ..., "chrom": ..., "start": ..., "end": ...}, ...]
    """
    text = Path(summary_path).read_text()
    regions = []

    # Find all repo paths mentioned
    repo_blocks = re.split(r"=== Test \(", text)[1:]   # one block per repo

    for block in repo_blocks:
        # extract repo path
        repo_match = re.match(r"(.+?)\) bin ranges by chromosome ===", block)
        if not repo_match:
            continue
        repo_path = repo_match.group(1).strip()

        # extract chromosome table rows
        lines = block.split("\n")
        for line in lines:
            # skip header lines
            if "bins_count" in line or "===" in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                chrom = parts[0]
                try:
                    min_start = int(parts[2])
                    max_end   = int(parts[3])
                    regions.append({
                        "repo": repo_path,
                        "chrom": chrom,
                        "start": min_start,
                        "end": max_end,
                    })
                except ValueError:
                    continue

    return regions

test_regions = parse_test_regions(SUMMARY_FILE)
if not test_regions:
    raise ValueError(f"Could not parse any test regions from {SUMMARY_FILE}")

# ============================================================
# Decide which repository to use
# ============================================================
# Detect all species present across all test regions
species_in_test = set()
for region in test_regions:
    sp = detect_species(region["chrom"])
    species_in_test.add(sp)

print(f"\nSpecies detected across all test regions: {species_in_test}")

# If more than one species is found → use concatenated momics
if len(species_in_test) > 1:
    DEFAULT_REPO = SPECIES_REPO["concatenated"]
    print(f"  → Multiple species detected, using concatenated repo: {DEFAULT_REPO}")
else:
    single_species = list(species_in_test)[0]
    DEFAULT_REPO = SPECIES_REPO.get(single_species, None)
    print(f"  → Single species detected ({single_species}), using repo: {DEFAULT_REPO}")


print(f"\nParsed {len(test_regions)} test region(s):")
for r in test_regions:
    print(f"  {r['repo']} | {r['chrom']}:{r['start']}-{r['end']}")

# ============================================================
# Load model
# ============================================================
print(f"\nLoading model from {MODEL_PATH}…")

custom_objects = {
    "Conv1DBlock":    nn_changed.Conv1DBlock,
    "DilatedConvBlock": nn_changed.DilatedConvBlock,
    "SoftGate":       nn_changed.SoftGate,
    "loss_mae_cor":   nn_changed.loss_mae_cor,
    "cor":            nn_changed.cor,
}

model = tf.keras.models.load_model(
    str(MODEL_PATH),
    custom_objects=custom_objects,
    compile=False,
)


# Ensure model is built to expose output shape
if model.output_shape is None:
    try:
        dummy = tf.zeros((1,) + model.input_shape[1:])
        _ = model(dummy)
    except Exception as e:
        print("Building model dynamically failed:", e)
    print("Model built to populate shape info.")

FEATURES_SIZE = model.input_shape[1]
# Handle sequential and multi-output models
out_shape = model.output_shape
if isinstance(out_shape, list):
    TARGET_SIZE = out_shape[0][1]
else:
    TARGET_SIZE = out_shape[1]
TARGETS = list(model.output_names)

print(f"✔ Model loaded")
print(f"  Input shape  : {model.input_shape}")
print(f"  Targets      : {TARGETS}")
print(f"  Feature size : {FEATURES_SIZE}")
print(f"  Target size  : {TARGET_SIZE}")


# ============================================================
# Always use the concatenated repo (handles multi-species)
# ============================================================
REPO_PATH = DEFAULT_REPO
if not Path(REPO_PATH).exists():
    raise FileNotFoundError(f"Repository not found: {REPO_PATH}")

print(f"\nUsing repository: {REPO_PATH}")
repo = Momics(REPO_PATH)
chrom_sizes = repo.chroms(as_dict=True)

# ============================================================
# Collect bins from ALL test regions into one combined set
# ============================================================
print(f"\nCollecting bins from all {len(test_regions)} test region(s)...")


all_bins_list = []

# --- Generate bins ONCE for the whole genome ---
print("Generating genome-wide bins (once)...")
all_repo_bins = repo.bins(
    width=FEATURES_SIZE,
    stride=args.stride,
    cut_last_bin_out=True
)
print(f"  Total genome bins available: {len(all_repo_bins)}")

for region in test_regions:
    CHROM  = region["chrom"]
    START  = region["start"]
    END    = region["end"]

    species = detect_species(CHROM)
    print(f"  {CHROM}:{START}-{END} → species: {species}")

    # --- Filter from pre-generated bins (no re-computation) ---
    region_bins = all_repo_bins[
        (all_repo_bins.Chromosome == CHROM)
        & (all_repo_bins.Start >= START - FEATURES_SIZE)
        & (all_repo_bins.End   <= END   + FEATURES_SIZE)
    ].copy()

    print(f"    → {len(region_bins)} bins found")

    if len(region_bins) > 0:
        all_bins_list.append(region_bins.df)


if not all_bins_list:
    raise ValueError("No bins found for any test region — check chromosome names and coordinates.")

# Combine all region bins into one DataFrame
import pyranges as pr
combined_bins_df = pd.concat(all_bins_list).reset_index(drop=True)
combined_bins = pr.PyRanges(combined_bins_df)
print(f"\nTotal combined bins: {len(combined_bins_df)}")

# ============================================================
# Streamer + TF dataset on combined bins
# ============================================================
print("\nCreating data streamer...")

streamer = MomicsStreamer(
    repo,
    combined_bins,
    features=["nucleotide"],
    batch_size=args.batch_size
)

ds = tf.data.Dataset.from_generator(
    lambda: streamer,
    output_signature={
        "nucleotide": tf.TensorSpec(
            shape=(None, FEATURES_SIZE, 4), dtype=tf.float32)
    },
)

# ============================================================
# Predict once on all bins
# ============================================================
print("\nRunning predictions...")
steps = int(np.ceil(len(combined_bins_df) / args.batch_size))
predictions = model.predict(ds, steps=steps, verbose=1)

# squeeze trailing dim if present
for tgt, arr in predictions.items():
    if arr.ndim == 3 and arr.shape[-1] == 1:
        predictions[tgt] = np.squeeze(arr, axis=-1)

print(f"✔ Predictions complete")
print(f"  Predicted tracks: {list(predictions.keys())}")

# ✅ CHECK 1 — are raw model predictions non-zero?
print("\n=== CHECK 1: Raw prediction value ranges ===")
for tgt, arr in predictions.items():
    print(f"  {tgt}: shape={arr.shape} | min={arr.min():.6f} | max={arr.max():.6f} | mean={arr.mean():.6f}")


# ============================================================
# Map predictions to genomic coordinates
# ============================================================
print("\nMapping predictions to genome coordinates...")

offset = (FEATURES_SIZE - TARGET_SIZE) // 2
bins_out = combined_bins_df.copy()
bins_out["Start"] = bins_out["Start"] + offset
bins_out["End"]   = bins_out["Start"] + TARGET_SIZE

# ✅ CHECK 2 — are bins_out coordinates and chromosomes correct?
print("\n=== CHECK 2: bins_out coordinates ===")
print(f"  Shape       : {bins_out.shape}")
print(f"  Chromosomes : {bins_out['Chromosome'].unique()}")
print(f"  Start range : {bins_out['Start'].min():,} → {bins_out['Start'].max():,}")
print(f"  End range   : {bins_out['End'].min():,} → {bins_out['End'].max():,}")
print(bins_out.head(3))

results = {}
for target in TARGETS:
    pred_key = find_prediction_key(target, predictions)
    if pred_key is None:
        print(f"  ⚠️ '{target}' not found in predictions "
              f"(available: {list(predictions.keys())}), skipping")
        continue
    print(f"  Matched '{target}' → '{pred_key}'")
    track_dict = {}
    for i, (chrom, start, end) in enumerate(zip(
            bins_out["Chromosome"], bins_out["Start"], bins_out["End"])):
        track_dict[f"{chrom}:{start}-{end}"] = predictions[pred_key][i]
    # use clean name (strip _softplus etc)
    clean_name = re.sub(r'_(softplus|head|output)$', '', target)
    results[clean_name] = track_dict


    # ✅ CHECK 3 — are track_dict keys and values correct?
    print(f"\n=== CHECK 3: track_dict sample for '{clean_name}' ===")
    sample_items = list(track_dict.items())[:3]
    for k, v in sample_items:
        arr_v = np.array(v)
        print(f"  key : '{k}'")
        print(f"  val : shape={arr_v.shape} | min={arr_v.min():.6f} | max={arr_v.max():.6f}")

# ============================================================
# Single output directory — one BigWig per track
# ============================================================
OUTPUT_DIR = Path(args.base_pred_dir) / args.run_name
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"\nSaving BigWigs to {OUTPUT_DIR}...")
original_dir = os.getcwd()
os.chdir(OUTPUT_DIR)

for target, track_dict in results.items():
    print(f"  Aggregating {target}...")
    mma.aggregate(
        {target: track_dict},
        chrom_sizes,
        type="mean",
        prefix=f"{target}",
    )

os.chdir(original_dir)
print(f"\n✅ All done. BigWigs saved to {OUTPUT_DIR}")

print(f"\n✅ All regions complete.")