#!/usr/bin/env python3
"""
Setup script for Momics-based training and inference.
Configures environment, TensorFlow settings, and imports
required tools and utilities.
"""

# =========================================================
# 1. Standard library imports
# =========================================================
import os
import sys
import argparse
import time
import csv
import gc
import json
import resource
from pathlib import Path
import psutil
import pandas as pd
import numpy as np
import pyranges as pr  # if not already imported
import pyBigWig

# =========================================================
# 2. Append Momics repository to the Python path
# =========================================================
sys.path.insert(0, "/pasteur/appa/scratch/adecugis/all_momics/momics/src")

# =========================================================
# 3. Resource limits (increase open files)
# =========================================================
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, 65536), hard))

# =========================================================
# 4. Environment variables (TileDB, TensorFlow, TMP)
# =========================================================
os.environ.update({
    "TILEDB_SM_CONSOLIDATION_MODE": "fragments",
    "TILEDB_SM_DEDUP_COORDS": "false",
    "TILEDB_VFS_FILE_MAX_PARALLEL_OPS": "1",
    "TILEDB_DISABLE_FILE_LOCKING": "1",
    "TILEDB_CONFIG": "",
    "TMPDIR": "/tmp",
    # TensorFlow GPU tuning
    "TF_GPU_THREAD_MODE": "gpu_private",
    "TF_GPU_THREAD_COUNT": "2",
    "TF_ENABLE_WINOGRAD_NONFUSED": "1",
    "TF_ENABLE_CUDNN_FRONTEND": "1",
})

# =========================================================
# 5. TensorFlow setup and session cleanup
# =========================================================
import tensorflow as tf
from tensorflow.keras import layers
from contextlib import redirect_stdout

# Clear previous TF session
try:
    import tensorflow.keras.backend as K
    K.clear_session()
except Exception:
    pass

gc.collect()

# Enable XLA JIT compilation
tf.config.optimizer.set_jit(True)

# Configure GPU memory growth
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    print(f"✓ Configured {len(gpus)} GPU(s)")

# =========================================================
# 6. Keras training callbacks
# =========================================================
from tensorflow.keras.callbacks import (
    CSVLogger,
    EarlyStopping,
    ModelCheckpoint,
    ReduceLROnPlateau,
)
import gc
import tensorflow.keras.backend as K
from tensorflow.keras.models import load_model
from momics.nn_changed_window import Conv1DBlock, DilatedConvBlock, loss_mae_cor  # import registers the class

# =========================================================
# 7. Mixed precision policy (for performance on GPU)
# =========================================================
from tensorflow.keras import mixed_precision
mixed_precision.set_global_policy('float32')
# =========================================================
# 8. Momics utilities and modules
# =========================================================

# from util functions
from momics.utils import (
    make_nan_replacer,
    prefilter_bins_for_nans_and_saves,
    ResourceLogger,
    ValidationSummaryCallback,
    EarlyStoppingAfterLRDrop,
)


import momics.nn_changed_window as nn_changed
from momics.momics import Momics
from momics.dataset import MomicsDataset


# =========================================================
# Argument parser placeholder 
# =========================================================

"""
takes in 3 argument
--momics _insert paths to momics repositories _
--tracks list_of_tracks_want_predict
--architecture _state your trunk of choice_

python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/momics_normalized_scale/cerevisiae.momics \
  --tracks ATAC STP1 \
  --architecture no_maxpool \
  --prefix_outdir june_12_STP1_ATAC \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/ \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/tracks_polymerase_normalized



python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/pombe.momics \
  --tracks ADN2 ATF1 BDF2 CUF1 FIL1 GAF1 H3K4ME3 H3K36ME3 HSR1 KLF1 MNASE PCR1 PHP3 PHP5 PHX1 POB3 REB1 RPB1 RSV2 SWR1 TBP1 \
  --architecture no_maxpool \
  --prefix_outdir better_bg_stdva6_june1 \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/Pombe_modulable \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe 

   
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/c_glabrata.momics \
  --tracks H3 MNASE_15 MNASE_25 MNASE_35 \
  --architecture no_maxpool \
  --prefix_outdir Manon_june2_cerevisiae_only \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/ \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \
  --test

  
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/cerevisiae.momics \
  --tracks ASH1 BRE1 BRE2 CBC2 CCL1 ABF1 H3K4ME3 H3K27AC H3K36ME3 MNASE RPB1 SCC1 STP1 CET1 CSN12 CTR9 CYC8 EAF7 ELF1 GLN3 HDA1 HMS2 HOM6 HPA3 SUA7 \
  --architecture no_maxpool \
  --prefix_outdir cerevesisiae_TEST_alpha30 \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/Cerevisiae_modulable \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input \
  --test

  python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/pombe.momics \
  --tracks BDF2 CUF1 GAF1 H3K4ME3 H3K36ME3 HSR1 KLF1 MNASE PCR1 PHP3 PHP5 PHX1 POB3 REB1 RPB1 RSV2 SWR1 TBP1 \
  --architecture no_maxpool_big \
  --prefix_outdir alpha20_27_may \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/Pombe_modulable \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \

python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/polymerase.momics \
  --tracks RPA14 RPA34 RPA43 RPA49 RPA135 RPA190 RPB2 RPB11 RPO21 RPC17 RPC25 RPC31 RPC34 RPC37 RPC53 RPC82 RPO31 \
  --architecture no_maxpool \
  --prefix_outdir polymerase_all_tracks_may26 \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/polymerase \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/tracks_polymerase_normalized

python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/polymerase.momics \
  --tracks RPA190 RPO21 RPO31 \
  --architecture no_maxpool \
  --prefix_outdir polymerase_3_tracks_may28_maecor_only \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/polymerase \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/tracks_polymerase_normalized


python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/cerevisiae.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/glabrata.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/pombe.momics \
  --tracks H3K4ME3 H3K36ME3 MNASE \
  --architecture no_maxpool \
  --prefix_outdir figures_methylation_multispecies \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/figures \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \
  
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/cerevisiae.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/glabrata.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/pombe.momics \
  --tracks H3K4ME3 \
  --architecture no_maxpool \
  --prefix_outdir figures_one_track_multispecies \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/figures \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \

"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--momics", nargs="+", type=str, required=True, help="paths to Momics repositories")
parser.add_argument("--tracks", nargs="+", type=str, required=True, help="list of track names to predict")
parser.add_argument("--architecture", type=str, default="default", help="Trunk architecture (e.g. default, kernel)")
parser.add_argument("--prefix_outdir", type=str, default="default", help="prefix_for_outdir")
parser.add_argument("--outdir_path", type=str, default="default", help="path for output directory")
parser.add_argument("--input_bw", type=str, default="/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe", help="path to BigWig files")
parser.add_argument("--test", action="store_true", help="Run 1 epoch of stage 1 only to test stage 2 transition")


args = parser.parse_args()

momics_paths_list = args.momics        # → list of strings
tracks_names_list = args.tracks        # → list of strings
trunk_name = args.architecture         # → single string
outdir_prefix = args.prefix_outdir            # → single string
outdir_path = args.outdir_path                      # → single string
bw_base_dir = Path(args.input_bw)
test_mode   = args.test         
print(f"test_mode = {test_mode}")   # ← should print True or False

logdir = Path(
    f"{outdir_path}_{outdir_prefix}_{trunk_name}"
)
logdir.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Repository + constants
# -----------------------------

#add other options
features_size = 8192 #20000  #8192 32_768
target_size = 512
stride = 48
batch_size = 20

#TF
alpha_input = 20
bg_input = 0.4

#methylation
alpha = 0.5 #alpha is cor
beta = 0.05


bg_method = "median_6stdv"  #if you change this name it will then reclculate even if you had saved the bed files before

#--------------------------
# Trunk registry
#--------------------

TRUNK_REGISTRY = {
    "no_maxpool": nn_changed.ChromNNBackboneUNet_TF_no_maxpool,
    "no_maxpool_big": nn_changed.ChromNNBackboneUNet_TF_no_maxpool_big
}


#-----------------
# Helper functions
#--------------

SPECIES_PATTERNS = {
    "CBS138":  "C_glabrata",
    "S288C":   "S_cerevisiae",
    "CDS216":  "S_eubayanus",
    "CBS432":  "S_paradoxus",
    "972h":    "S_pombe",
    "CU329":   "S_pombe",
}

def detect_species(chrom_name: str) -> str:
    """
    Map a chromosome name to a species ID using SPECIES_PATTERNS.
    Returns 'unknown' if no match found.
    """
    for pattern, species in SPECIES_PATTERNS.items():
        if pattern in chrom_name:
            return species
        
    return "unknown"


def save_bg_clip_passing_bins(bins_df, bg_clip_values, bw_base_dir, output_dir, bg_method):
    output_dir  = Path(output_dir)
    bw_base_dir = Path(bw_base_dir)

    bed_output_dir = bw_base_dir / "bed_files_pass_bg_clip" / bg_method
    bed_output_dir.mkdir(parents=True, exist_ok=True)

    bg_clip_txt = bed_output_dir / "bg_clip_values.txt"
    with open(bg_clip_txt, "w") as f:
        f.write("track\tbg_clip\n")
        for track, bg_clip in bg_clip_values.items():
            f.write(f"{track}\t{bg_clip:.6f}\n")
    print(f"  ✔ Saved bg_clip values → {bg_clip_txt}")

    for track, bg_clip in bg_clip_values.items():

        # ── ADD HERE: skip if already exists ─────────────────────  ← NEW
        bed_path = bed_output_dir / f"{track}_bg_clip_passing_bins.bed" 
        
        if bed_path.exists():                                            
            print(f"  ⏭️  {track}: BED already exists, skipping")        
            continue                                                     

        bw_file = bw_base_dir / f"{track}.bw"
        if not bw_file.exists():
            print(f"  ⚠️  BigWig not found for {track}, skipping BED")
            continue

        passing = []
        bw = pyBigWig.open(str(bw_file))

        for _, row in bins_df.iterrows():
            chrom = row["Chromosome"]
            start = int(row["Start"])
            end   = int(row["End"])
            try:
                arr = np.array(bw.values(chrom, start, end, numpy=True))
                arr = arr[np.isfinite(arr)]
                if arr.size > 0 and arr.max() > bg_clip:
                    passing.append((chrom, start, end))
            except Exception:
                continue

        bw.close()

        # ── bed_path already defined above, just write ────────────  ← MOVED
        with open(bed_path, "w") as f:
            for chrom, start, end in passing:
                f.write(f"{chrom}\t{start}\t{end}\n")

        n_total = len(bins_df)
        n_pass  = len(passing)
        pct     = 100 * n_pass / n_total if n_total > 0 else 0
        print(f"  {track}: {n_pass}/{n_total} bins ({pct:.1f}%) "
              f"pass bg_clip={bg_clip:.4f} → {bed_path}")
        
def save_run_config(model, logdir, args, bg_clip_values):
    """
    Save a JSON file documenting:
    - CLI arguments
    - All training constants
    - Per-track bg_clip values
    - Per-layer architecture details (type, kernel size, filters, params)
    - Full model summary as text
    """

    # ── constants ──────────────────────────────────────────────────
    constants = {
        "features_size": features_size,
        "target_size":   target_size,
        "stride":        stride,
        "batch_size":    batch_size,
        "alpha_input":   alpha_input,
        "bg_input":      bg_input,
        "alpha":         alpha,
        "beta":          beta,
    }

    # ── CLI / run info ──────────────────────────────────────────────
    run_info = {
        "trunk_name":    args.architecture,
        "prefix_outdir": args.prefix_outdir,
        "outdir_path":   args.outdir_path,
        "momics_paths":  args.momics,
        "tracks":        args.tracks,
    }

    # ── per-layer architecture details ─────────────────────────────
    arch_layers = []
    for layer in model.layers:
        try:
            cfg = layer.get_config()
        except Exception:
            cfg = {}

        # pull out the most useful fields cleanly
        arch_layers.append({
            "name":             layer.name,
            "type":             type(layer).__name__,
            "trainable_params": int(layer.count_params()),
            "trainable":        bool(layer.trainable),
            # key architectural numbers (present only for Conv layers)
            "filters":          cfg.get("filters"),
            "kernel_size":      cfg.get("kernel_size"),
            "strides":          cfg.get("strides"),
            "dilation_rate":    cfg.get("dilation_rate"),
            "padding":          cfg.get("padding"),
            "activation":       cfg.get("activation"),
            "full_config":      cfg,
        })

    # ── model summary as list of strings ───────────────────────────
    summary_lines = []
    model.summary(print_fn=lambda line: summary_lines.append(line))

    # ── assemble ───────────────────────────────────────────────────
    run_config = {
        "run_info":         run_info,
        "constants":        constants,
        "bg_clip_values":   {k: float(v) for k, v in bg_clip_values.items()},
        "total_params":     int(model.count_params()),
        "architecture":     arch_layers,
        "model_summary":    summary_lines,
    }

    json_path = Path(logdir) / "run_config.json"
    with open(json_path, "w") as f:
        json.dump(run_config, f, indent=2, default=str)

    print(f"✔ Saved run config → {json_path}")

def compute_bg_clip_from_bw(bw_path):
    """
    Compute bg_clip = median(signal) + 2 * std(signal)
    """
    bw_path = str(bw_path)

    if not os.path.exists(bw_path):
        print(f"  ⚠️ File not found: {bw_path}")
        return 0.4

    bw = pyBigWig.open(bw_path)
    vals = []

    for chrom in bw.chroms().keys():
        arr = np.array(bw.values(chrom, 0, bw.chroms()[chrom], numpy=True))
        arr = arr[np.isfinite(arr)]
        if arr.size > 0:
            vals.append(arr)

    bw.close()
    if not vals:
        return 0.0

    vals = np.concatenate(vals)
    bg_clip = np.median(vals) + 6 * np.std(vals) #change number of stdva for clip

    #if (np.median(vals) + 7 * np.std(vals)) > 0.4:
    #    bg_clip = 0.4
    #else:
    #    bg_clip = np.median(vals) + 7 * np.std(vals)

    return bg_clip


def split_training_data_stratified(
    momics_paths_list,
    features_size,
    target_size,
    stride,
    batch_size,
    split_ratio=(0.8, 0.1, 0.1)
):
    """
    For each Momics repository, split bins into train/val/test
    ensuring each species contributes proportionally.
    Bins are sorted by position so each split is a continuous
    genomic region (not random).
    """
    split_summary = {}
    min_train_bins = float('inf')
    corresponding_val_bins = None

    for repo_path_sub in momics_paths_list:
        print(f"\n{'='*60}")
        print(f"Processing repository: {repo_path_sub}")
        print(f"{'='*60}")

        repo = Momics(repo_path_sub)
        bins = repo.bins(
            width=features_size,
            stride=stride,
            cut_last_bin_out=True
        )

        bins_df = bins.df.copy()

        # --- detect species per bin ---
        bins_df["species"] = bins_df["Chromosome"].apply(detect_species)
        species_found = bins_df["species"].unique()
        print(f"Species found in repo: {list(species_found)}")

        train_parts, val_parts, test_parts = [], [], []

        for sp in species_found:
            sp_df = bins_df[bins_df["species"] == sp].copy()
            n = len(sp_df)
            n_train = int(n * split_ratio[0])
            n_val   = int(n * split_ratio[1])
            n_test  = n - n_train - n_val

            # sort by position → continuous blocks
            sp_df = sp_df.sort_values(["Chromosome", "Start"]).reset_index(drop=True)

            train_parts.append(sp_df.iloc[:n_train])
            val_parts.append(sp_df.iloc[n_train: n_train + n_val])
            test_parts.append(sp_df.iloc[n_train + n_val:])

            # report the actual genomic coordinates of each split
            train_block = sp_df.iloc[:n_train]
            val_block   = sp_df.iloc[n_train: n_train + n_val]
            test_block  = sp_df.iloc[n_train + n_val:]

            print(f"\n  {sp}: total={n} | train={n_train} | val={n_val} | test={n_test}")
            if len(train_block) > 0:
                print(f"    Train : {train_block.iloc[0]['Chromosome']} "
                      f"{train_block.iloc[0]['Start']:,} → "
                      f"{train_block.iloc[-1]['End']:,}")
            if len(val_block) > 0:
                print(f"    Val   : {val_block.iloc[0]['Chromosome']} "
                      f"{val_block.iloc[0]['Start']:,} → "
                      f"{val_block.iloc[-1]['End']:,}")
            if len(test_block) > 0:
                print(f"    Test  : {test_block.iloc[0]['Chromosome']} "
                      f"{test_block.iloc[0]['Start']:,} → "
                      f"{test_block.iloc[-1]['End']:,}")

        # combine across species
        train_df = pd.concat(train_parts).reset_index(drop=True)
        val_df   = pd.concat(val_parts).reset_index(drop=True)
        test_df  = pd.concat(test_parts).reset_index(drop=True)

        n_train = len(train_df)
        n_val   = len(val_df)
        n_test  = len(test_df)

        print(f"\nTotal → train={n_train} | val={n_val} | test={n_test}")

        split_summary[repo_path_sub] = {
            "train_df": train_df,
            "val_df":   val_df,
            "test_df":  test_df,
            "n_train":  n_train,
            "n_val":    n_val,
            "n_test":   n_test,
        }

        if n_train < min_train_bins:
            min_train_bins = n_train
            corresponding_val_bins = n_val

    print("\n" + "="*60)
    print(f"Minimum training bins across repos: {min_train_bins}")
    print(f"Corresponding val bins: {corresponding_val_bins}")
    print("="*60)

    return split_summary, min_train_bins, corresponding_val_bins


def summarize_ranges(df, name="Dataset"):
    summary = (
        df.groupby("Chromosome")
        .agg(
            bins_count=("Chromosome", "count"),
            min_start=("Start", "min"),
            max_end=("End", "max")
        )
        .sort_index()
    )
    print(f"\n=== {name} bin ranges by chromosome ===")
    print(summary)

def build_multi_repo_datasets(
    repo_paths,
    features, #"features": ["nucleotide"] as default
    targets,
    features_size,
    target_size,
    stride,
    batch_size,
    min_train_bins,             #  mandatory: number of training bins to keep per repo
    corresponding_val_bins,    #  mandatory: number of validation bins to keep per repo
    split_summary, 
    outdir_path,
):
    """
    Build balanced multi‑species Momics datasets using a fixed number of bins
    (in_train_bins and corresponding_val_bins).  The remaining bins per
    species are returned as the coordinates for its test region.
    """

    #set features default
    if features is None:
        features = ["nucleotide"]
    if targets is None:
        targets = []

    train_datasets = []
    val_datasets = []
    test_regions  = {}
    total_train_bins = 0

    for repo_path_sub in repo_paths:
        print(f"\n{'='*60}")
        print(f"Processing repository: {repo_path_sub}")
        print(f"{'='*60}")

        repo = Momics(repo_path_sub)

        # Use pre-computed stratified splits
        repo_split = split_summary[repo_path_sub]
        train_bins_df = repo_split["train_df"].drop(columns=["species"], errors="ignore")
        val_bins_df   = repo_split["val_df"].drop(columns=["species"], errors="ignore")
        test_bins_df  = repo_split["test_df"].drop(columns=["species"], errors="ignore")

        # Balance: cap training and val to minimum across species
        train_bins_df = train_bins_df.iloc[:min_train_bins]
        val_bins_df   = val_bins_df.iloc[:corresponding_val_bins]

        print(f"Stratified split → train: {len(train_bins_df)} | val: {len(val_bins_df)} | test: {len(test_bins_df)}")

        summarize_ranges(train_bins_df, "Train")
        summarize_ranges(val_bins_df,   "Validation")
        summarize_ranges(test_bins_df,  "Test")

        # Convert to PyRanges
        train_bins = pr.PyRanges(train_bins_df[["Chromosome","Start","End"]])
        val_bins   = pr.PyRanges(val_bins_df[["Chromosome","Start","End"]])
        test_bins  = pr.PyRanges(test_bins_df[["Chromosome","Start","End"]])

        # NaN filtering
        print("Filtering training bins for NaNs…")
        train_bins = prefilter_bins_for_nans_and_saves(
            repo=repo, bins=train_bins, targets=targets,
            max_empty_bins=20, target_size=target_size)

        print("Filtering validation bins for NaNs…")
        val_bins = prefilter_bins_for_nans_and_saves(
            repo=repo, bins=val_bins, targets=targets,
            max_empty_bins=20, target_size=target_size)

        # Store test bins per species for prediction later
        test_regions[repo_path_sub] = test_bins[["Chromosome","Start","End"]]

        # Build TF datasets
        train_ds = MomicsDataset(
            repo, train_bins, features, targets,
            target_size=target_size, batch_size=batch_size)
        val_ds = MomicsDataset(
            repo, val_bins, features, targets,
            target_size=target_size, batch_size=batch_size)

        train_datasets.append(train_ds)
        val_datasets.append(val_ds)
        total_train_bins += len(train_bins)

    # Combine, preprocess, prefetch
    train_dataset = tf.data.Dataset.sample_from_datasets(train_datasets)
    val_dataset   = tf.data.Dataset.sample_from_datasets(val_datasets)
    nan_replacer  = make_nan_replacer(targets)

    train_dataset = (train_dataset.unbatch()
                     .map(nan_replacer, num_parallel_calls=tf.data.AUTOTUNE)
                     .batch(batch_size, drop_remainder=True)
                     .repeat()
                     .prefetch(tf.data.AUTOTUNE))

    val_dataset = (val_dataset.unbatch()
                   .map(nan_replacer, num_parallel_calls=tf.data.AUTOTUNE)
                   .batch(batch_size, drop_remainder=True)
                   .prefetch(tf.data.AUTOTUNE))

    # Save test region summary
    log_path = os.path.join(logdir, "test_region_summary.txt")
    with open(log_path, "w") as f:
        print("\n=== Test regions (remaining bins) ===", file=f)
        for repo_p, region_pr in test_regions.items():
            print(f"{repo_p}: {len(region_pr)} test bins", file=f)
            with redirect_stdout(f):
                summarize_ranges(region_pr.df, f"Test ({repo_p})")
    print(f"Saved summary to: {log_path}")

    return train_dataset, val_dataset, total_train_bins, test_regions

# -----------------------------
# Training function
# -----------------------------

def train_model(tracks_list: list[str], trunk_name: str):
    
    # ── bg_clip values ───────────────────────────────────────────
    bg_clip_values = {}
    for track in tracks_names_list:
        bw_file = bw_base_dir / f"{track}.bw"
        if bw_file.exists():
            clip_val = compute_bg_clip_from_bw(bw_file)
            bg_clip_values[track] = clip_val
            print(f"{track}: bg_clip = {clip_val:.4f}")
        else:
            print(f"⚠️ bigWig file for {track} not found, using default 0.4")
            bg_clip_values[track] = 0.4

    # ── ✅ NEW: BED files for all bins passing bg_clip per track ──
    print("\nGenerating bg_clip BED files …")
    for repo_path_sub in momics_paths_list:
        repo      = Momics(repo_path_sub)
        all_bins  = repo.bins(width=features_size, stride=stride, cut_last_bin_out=True)
        save_bg_clip_passing_bins(
            bins_df       = all_bins.df,
            bg_clip_values= bg_clip_values,
            bw_base_dir   = bw_base_dir,
            output_dir    = logdir / "bg_clip_beds",
            bg_method  = bg_method,
        )

    HISTONE_MARKS = []
    TF = []



    #split data
    split_summary, min_train_bins, corresponding_val_bins =  split_training_data_stratified(
        momics_paths_list=momics_paths_list,
        features_size=features_size,
        target_size=target_size,
        stride=stride,
        batch_size=batch_size,
        split_ratio=(0.8, 0.1, 0.1))

    #sort targets in 2 groups
    for t in tracks_list:
        t_upper = t.upper()
        if t_upper.startswith("MNASE") or t_upper.startswith("H3") or t_upper in {"RPO31", "RPA190", "RPO21", "ATAC", "STP1"}:
            HISTONE_MARKS.append(t)
        else:
            TF.append(t)

    def make_loss_mae_cor_weighted(alpha=alpha, beta=beta):  #alpha is cor
        def loss_fn(y_true, y_pred):
            return nn_changed.loss_mae_cor_weighted(y_true, y_pred, alpha=alpha, beta=beta)
        return loss_fn

    def make_mae_loss(alpha, bg_clip):
        def loss_fn(y_true, y_pred):
            return nn_changed.weighted_mae(y_true, y_pred, alpha=alpha, bg_clip=bg_clip)
            #return nn_changed.loss_mae_cor_weighted(y_true, y_pred, alpha=alpha, beta=beta)
        return loss_fn

    def get_loss_for_target(t):
        if t in HISTONE_MARKS:
            fn = make_loss_mae_cor_weighted(alpha=alpha, beta=beta)
            fn.__name__ = f"loss_mae_cor_weighted(alpha={alpha}, beta={beta})"
            return fn
        elif t in TF:
            clip_val = bg_clip_values.get(t, bg_input)
            fn = make_mae_loss(alpha=alpha_input, bg_clip=clip_val)
            fn.__name__ = f"loss_mae(alpha={alpha_input}, clip={clip_val:.4f})"
            return fn
        else:
            raise ValueError(f"Target '{t}' not found in any group. Please add it.")
        

    # Build per-head loss dict
    loss_dict = {}

    print("Loss function assignment per target:")

    for t in tracks_list:
        loss_fn = get_loss_for_target(t)
        loss_dict[t] = loss_fn
        print(f"  {t}: {getattr(loss_fn, '__name__', str(loss_fn))}")

    # ---------------------------------------------------------
    # Build datasets
    # ---------------------------------------------------------

    train_dataset, val_dataset, bins_train_number, test_regions = build_multi_repo_datasets(
        repo_paths=momics_paths_list,
        features=None,              # replace with actual feature list when available
        targets=tracks_list,  # from argparse or define above
        features_size=features_size,
        target_size=target_size,
        stride=stride,
        batch_size=batch_size,
        min_train_bins=min_train_bins,
        corresponding_val_bins=corresponding_val_bins,
        split_summary=split_summary,
        outdir_path=logdir,  # pass logdir for any saving needs within the function
    )

    # ---------------------------------------------------------
    # Choose trunk architectures --> head is included
    # ---------------------------------------------------------
    trunk_cls = TRUNK_REGISTRY[trunk_name]

    # ---------------------------------------------------------
    # Build model 1
    # ---------------------------------------------------------
    inputs = {"nucleotide": layers.Input(shape=(features_size, 4), name="nucleotide"),}
    outputs = {t: layers.Reshape((target_size,)) for t in tracks_list}

    model = nn_changed.ChromNNBase_modulable(
        inputs=inputs,
        outputs=outputs,
        trunk_cls=lambda inputs: trunk_cls(inputs, features_size=features_size, target_size=target_size),
        name=f"ChromNN_{trunk_name}",
    ).model

    print("Model outputs:", list(model.output_names))
    print("Loss dict keys:", list(loss_dict.keys()))  # if you built one per target

    metrics_list = [
        nn_changed.spearman_cor,
        nn_changed.pearson_cor,
        nn_changed.weighted_mae,
        nn_changed.loss_mae_cor,
    ]

    for k, v in loss_dict.items():
        print(k, v)

    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4, global_clipnorm=1.0),
        loss=loss_dict,
        metrics={t: metrics_list for t in tracks_list},
    )

    model.summary()

    # ── ✅ NEW: JSON run config (call right after model.summary()) ─
    save_run_config(
        model          = model,
        logdir         = logdir,
        args           = args,
        bg_clip_values = bg_clip_values,
    )

    validation_summary_cb = ValidationSummaryCallback(
        training_csv_path=str(logdir / "training.csv"),
        out_dir=str(logdir),
        tracks=tracks_list
    )   

    callbacks_stage1 = [
        CSVLogger(logdir / "training.csv", append=True),
        ResourceLogger(logdir / "resources.csv"),
        ReduceLROnPlateau(
            monitor="val_loss",
            factor=0.5,
            patience=2,
            min_delta=1e-4,
            min_lr=1e-6,
            verbose=1,
        ),
        EarlyStoppingAfterLRDrop(    # replaces EarlyStopping
            monitor="val_loss",
            patience=3,
            min_delta=1e-4,
            min_lr_factor=0.5,       # stop only after LR ≤ initial_lr * 0.5
            restore_best_weights=True,
        ),
        validation_summary_cb,
    ]
    

    #---------------------------------------------------------
    # Steps (override in test mode)
    # ---------------------------------------------------------

    if test_mode:
        steps_per_epoch = 2
        val_steps       = 2
        epochs_stage1   = 1          # just 1 epoch so we drop straight to stage 2
        print("\n TEST MODE: 1 epoch, 2 steps — checking stage 2 transition\n")
    else:
        steps_per_epoch = bins_train_number // batch_size
        val_steps       = (bins_train_number // batch_size) // 10
        epochs_stage1   = 150

    print(f"\nStarting stage 1: training trunk and head together...\n")
    print("steps per epoch:", steps_per_epoch)
    print("bins train number", bins_train_number)

    (logdir / "checkpoints_stage1").mkdir(parents=True, exist_ok=True)        # stage 1
    checkpoint_cb_all = ModelCheckpoint(
        filepath=logdir / "checkpoints_stage1/epoch_{epoch:03d}.keras",
        save_weights_only=False,
        save_freq="epoch"
    )

    checkpoint_cb_best = ModelCheckpoint(
        filepath=logdir / "trunk_model.keras",
        save_weights_only=False,
        save_best_only=True,
        monitor="val_loss",
        mode="min"
    )

    model.fit(
        train_dataset,
        validation_data=val_dataset,
        epochs=epochs_stage1,        # ← was hardcoded 150
        steps_per_epoch=steps_per_epoch,
        validation_steps=val_steps,
        callbacks=[checkpoint_cb_all, checkpoint_cb_best, *callbacks_stage1],
        verbose=1,
    )
    
    model.save(logdir / "trunk_model.keras")
    


    #---------------------------------------------------------
    # Clear memoory
    #---------------------------------------------------------

    # Clear current TensorFlow session
    K.clear_session()
    gc.collect()
    
    try:
        tf.config.experimental.reset_memory_stats('GPU:0')
    except Exception as e:
        print(f"  ⚠️ Could not reset memory stats: {e}")


    print("Cleared TensorFlow and Python memory, ready to reload stage 2 model.")

    validation_summary_cb = ValidationSummaryCallback(   # ← recreate here
        training_csv_path=str(logdir / "training.csv"),
        out_dir=str(logdir),
        tracks=tracks_list
    )

    # ---------------------------------------------------------
    # Build model 2
    # ---------------------------------------------------------

    model = load_model(
        logdir / "trunk_model.keras",
        custom_objects={
            "Conv1DBlock":           nn_changed.Conv1DBlock,
            "DilatedConvBlock":      nn_changed.DilatedConvBlock,
            "loss_mae_cor":          nn_changed.loss_mae_cor,
            "loss_mae_cor_weighted": nn_changed.loss_mae_cor_weighted,
            "weighted_mae":          nn_changed.weighted_mae,
            "spearman_cor":          nn_changed.spearman_cor,
            "pearson_cor":           nn_changed.pearson_cor,
        },
        compile=False,    # ← skip loss deserialization entirely
    )

    # Then manually recompile with the loss_dict you already built above
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-5, global_clipnorm=1.0),
        loss=loss_dict,       # ← reuse the dict built earlier in train_model
        metrics={t: metrics_list for t in tracks_list},
    )

    #stage 2 - train nothing head and recompile separatly after freezing
    for layer in model.layers:
        if "head" not in layer.name:
            layer.trainable = False
        else:
            layer.trainable = True

    #smaller learning rate for fine-tuning head with frozen trunk
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-5, global_clipnorm=1.0),
        loss=loss_dict,
        metrics={t: metrics_list for t in tracks_list},
    )

    model.summary()

    #second set of call callbacks with different early stopping settings

    callbacks_stage2 = [
        CSVLogger(logdir / "training.csv", append=True),
        ResourceLogger(logdir / "resources.csv"),
        EarlyStopping(
            monitor="val_loss",
            patience=3,
            min_delta=1e-5,
            restore_best_weights=True,
        ),
        ReduceLROnPlateau(
            monitor="val_loss",
            factor=0.5,
            patience=6,
            min_lr=1e-5,
            verbose=1,
        ),
        validation_summary_cb,  
    ]

    print("\nStarting stage 2: training head with frozen trunk...\n")

    (logdir / "checkpoints_stage2").mkdir(parents=True, exist_ok=True)
    checkpoint_cb_all = ModelCheckpoint(
        filepath=logdir / "checkpoints_stage2/epoch_{epoch:03d}.keras",
        save_weights_only=False,
        save_freq="epoch"
    )

    checkpoint_cb_best = ModelCheckpoint(
        filepath=logdir / "trunk_model_stage2.keras",
        save_weights_only=False,
        save_best_only=True,
        monitor="val_loss",
        mode="min"
    )

    if test_mode:
        epochs_stage2 = 1
        steps_per_epoch = 2
        val_steps       = 2
        print("\n TEST MODE: 1 epoch for stage 2\n")
    else:
        epochs_stage2 = 150
        steps_per_epoch = bins_train_number // batch_size
        val_steps       = (bins_train_number // batch_size) // 10

    model.fit(
        train_dataset,
        validation_data=val_dataset,
        epochs=epochs_stage2,
        steps_per_epoch=steps_per_epoch,
        validation_steps=val_steps,
        callbacks=[checkpoint_cb_all, checkpoint_cb_best, *callbacks_stage2],
        verbose=1,
    )
    
    model.save(logdir / "final_head_model.keras")

# =========================================================
# Main block
# =========================================================
if __name__ == "__main__":
    train_model(tracks_names_list, trunk_name)
