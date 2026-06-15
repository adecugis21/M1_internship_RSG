#!/bin/bash

#SBATCH --job-name=2s1t
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=20:00:00
#SBATCH --output=/pasteur/appa/scratch/adecugis/all_momics/momics/logs/chromnn_%j.out
#SBATCH --error=/pasteur/appa/scratch/adecugis/all_momics/momics/logs/chromnn_%j.err

set -e
set -o pipefail

source ~/.bashrc
micromamba activate momics_fresh

ulimit -n 65536
echo "File limit: $(ulimit -n)"

export TILEDB_DISABLE_FILE_LOCKING=1
export TILEDB_VFS_FILE_MAX_PARALLEL_OPS=1
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH


#Figure 2
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/cerevisiae.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/pombe.momics \
  --tracks H3K4ME3  \
  --architecture no_maxpool \
  --prefix_outdir june4_2_species_1tracks \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/figures \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \
  
#figure 1b TF
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_and_databases/c_glabrata.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_and_databases/S_cerevisiae.momics \
  --tracks GAF1 H3K4ME3 H3K36ME3 MNASE \
  --architecture no_maxpool \
  --prefix_outdir compare_to_betterbg_june1_1TF_training_fig1b_june3 \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/Pombe_modulable \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe 


python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/momics_objects/blasting_TF/pombe.momics \
  --tracks CUF1 GAF1 H3K4ME3 H3K36ME3 KLF1 MNASE \
  --architecture no_maxpool \
  --prefix_outdir better_bg_stdva6_june1 \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/Pombe_modulable \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe 

#figure 2b MNASE manon +dtaabase
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_and_databases/c_glabrata.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_and_databases/S_cerevisiae.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_and_databases/S_pombe.momics \
  --tracks MNASE MNASE_15 MNASE_25 MNASE_35 \
  --architecture no_maxpool \
  --prefix_outdir 4MNASE_june5_test \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/Pombe_modulable \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \
  --test


#figure 3a
python /pasteur/helix/scratch/adecugis/all_momics/momics/scripts/modulable_script/V4_normal_AND_concatenated.py \
  --momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/cerevisiae.momics /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/pombe.momics \
  --tracks H3K4ME3  \
  --architecture no_maxpool \
  --prefix_outdir june4_2_species_1tracks \
  --outdir_path /pasteur/appa/scratch/adecugis/all_momics/momics/output/figures \
  --input_bw /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/blasting_TF_input_pombe \
  
