#!/bin/bash

#SBATCH --job-name=tinymap
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=1:00:00

bowtie2-build ./S288c.fa ./S288c
bowtie2-build ./CBS138.fa ./CBS138
/pasteur/appa/scratch/adecugis/handle_files/inputs/inputs_actually_used_final/LM184_nxq_R1.fq.gz
#with output
srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample SRR24512033 \
   --input SRR24512032 \
   --genome /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_32_33


srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample SRR8719730 \
   --input SRR8719728 \
   --genome /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/SRR8719730

srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample SRR8719729 \
   --input SRR8719727 \
   --genome /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/SRR8719729





   --calibration /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138 \


#----------------5 yeasts------------------

cd /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/Data/

srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample MP299_nxq\
   --input MP302_nxq  \
   --genome /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/merged_yeasts \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_5_yeasts_H3


srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample MP335_nxq\
   --input MP302_nxq  \
   --genome /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/merged_yeasts \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_5_yeasts_MNASE_15min

srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample MP336_nxq\
   --input MP302_nxq  \
   --genome /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/merged_yeasts \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_5_yeasts_MNASE_25min

srun /pasteur/appa/homes/adecugis/bin/tinyMapper/tinyMapper.sh \
   --mode ChIP \
   --sample MP337_nxq\
   --input MP302_nxq  \
   --genome /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/merged_yeasts \
   --threads 8 \
   --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_5_yeasts_MNASE_35min


#----------------------------
cd /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs

bowtie2 -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138 \
  -1 SRR24895598_1.fastq \
  -2 SRR24895598_2.fastq \
  -S SRR24895598.sam \
  -p 8

prefetch SRR24512032 SRR24512033
fastq-dump --gzip --split-3 SRR24512032
fastq-dump --gzip --split-3 SRR24512033

---------------------------------------
module load bowtie2 samtools

bowtie2 -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference//pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/merged_yeasts.fa \
  -1 /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/Data/MP299_nxq_R1.fq.gz \
  -2 /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/Data/MP299_nxq_R2.fq.gz \
  | samtools sort -o MP299_nxq.bam -
samtools index MP299_nxq.bam

bowtie2 -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138 \
  -1 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SCC1_lab/CH227_S16_R1.fq.gz \
  -2 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SCC1_lab/CH227_S16_R2.fq.gz \
  | samtools sort -o CH227_S16.bam -
samtools index CH227_S16.bam

bowtie2 -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  -1 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR24512032_R1.fastq.gz \
  -2 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR24512032_R2.fastq.gz \
  | samtools sort -o SRR24512032.bam -
samtools index SRR24512032.bam

bowtie2 -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  -1 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR8719729_R1.fastq.gz \
  -2 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR8719729_R2.fastq.gz \
  | samtools sort -o SRR8719729.bam -
samtools index SRR8719729.bam

bowtie2 -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  -1 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR24512033_R1.fastq.gz \
  -2 /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR24512033_R2.fastq.gz \
  | samtools sort -o SRR24512033.bam -
samtools index SRR24512033.bam




#if one stranded --> make bw and divide
samtools index SRR948227_IP.sorted.bam
samtools index SRR948234_input.sorted.bam

#making bw from spiked bam
bamCoverage   --bam LM188_nxq^unmapped_S288c^mapped_CBS138^filtered^8TZSMB.bam   --outFileName LM184_nxq_CBS138_spikein.bw   --binSize 10   --normalizeUsing CPM   --numberOfProcessors 8

#get SRR file
fasterq-dump --split-files SRR8719727.sra
fasterq-dump --split-files SRR8719729
25(input) and 28 (IP)

gzip SRR948234.fastq

gzip *.fastq

for f in *.fastq; do
  # compress the file
  gzip "$f"
done

# then rename _1/_2 → _R1/_R2
for f in *_1.fastq.gz; do
  mv "$f" "${f/_1.fastq.gz/_R1.fastq.gz}"
done

for f in *_2.fastq.gz; do
  mv "$f" "${f/_2.fastq.gz/_R2.fastq.gz}"
done

#make ip vs input file bw
cd /pasteur/appa/homes/adecugis/handle_files/output_H3K36me3Glabrata/bam/spikein


/pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR8719727.fasta.gz
#Run Mnase
srun tinyMapper.sh \
  --mode ChIP \
  --sample /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR8719727.fasta.gz \
  --genome /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  --threads 8 \
  --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_SRR8719727_SCC1

srun tinyMapper.sh \
  --mode MNase \
  --sample /pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/Spombe_MNase \
  --genome /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  --threads 8 \
  --output /pasteur/appa/scratch/adecugis/handle_files/outputs/output_pombe_MNase



# ============================================
# 1. Process IP (SRR26031578)
# ============================================
cd SRR26031578

samtools view -bS SRR26031578^mapped_CBS138^P6DEZN.sam | \
  samtools sort -@ 8 -o SRR26031578_CBS138.sorted.bam

samtools index SRR26031578_CBS138.sorted.bam

bamCoverage \
  --bam SRR26031578_CBS138.sorted.bam \
  --outFileName SRR26031578_CBS138_spikein.bw \
  --binSize 10 \
  --normalizeUsing CPM \
  -p 8

cd ..

# ============================================
# 2. Process Input (SRR26031012)
# ============================================
cd SRR26031012

samtools view -bS SRR26031012^mapped_CBS138^P6DEZN.sam | \
  samtools sort -@ 8 -o SRR26031012_CBS138.sorted.bam

samtools index SRR26031012_CBS138.sorted.bam

bamCoverage \
  --bam SRR26031012_CBS138.sorted.bam \
  --outFileName SRR26031012_CBS138_spikein.bw \
  --binSize 10 \
  --normalizeUsing CPM \
  -p 8

cd ..

# ============================================
# 3. Calculate IP/Input Ratio (QC)
# ============================================
bigwigCompare \
  --bigwig1 SRR26031578/SRR26031578_CBS138_spikein.bw \
  --bigwig2 SRR26031012/SRR26031012_CBS138_spikein.bw \
  --operation log2 \
  --outFileName SRR26031578_over_SRR26031012_CBS138_log2ratio.bw \
  --binSize 10 \
  -p 8

echo "✔ All spike-in tracks created"
ls -lh */*.bw *.bw

#---------------------------------------------------------------------------


cd /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects

momics create S_cerevisiae.momics #
momics ingest chroms -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/chrom_sizes_renamed/S288c.chrom.sizes -g merged_yeasts S_cerevisiae.momics
momics ingest seq -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/ref_genomes_renamed/S288c.fa S_cerevisiae.momics

momics ingest tracks -f H3=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_cerevisiae/H3_concatenated_norm0999.bw S_cerevisiae.momics
momics ingest tracks -f MNASE_15=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_cerevisiae/MNASE_15min_concatenated_norm0999.bw S_cerevisiae.momics
momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_cerevisiae/MNASE_25min_concatenated_norm0999.bw S_cerevisiae.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_cerevisiae/MNASE_35min_concatenated_norm0999.bw S_cerevisiae.momics

momics create c_glabrata.momics #
momics ingest chroms -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/chrom_sizes_renamed/CBS138.chrom.sizes -g merged_yeasts c_glabrata.momics
momics ingest seq -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/ref_genomes_renamed/CBS138.fa c_glabrata.momics

momics ingest tracks -f H3=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/C_glabrata/H3_concatenated_norm0999.bw c_glabrata.momics
momics ingest tracks -f MNASE_15=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/C_glabrata/MNASE_15min_concatenated_norm0999.bw c_glabrata.momics
momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/C_glabrata/MNASE_25min_concatenated_norm0999.bw c_glabrata.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/C_glabrata/MNASE_35min_concatenated_norm0999.bw c_glabrata.momics

momics create S_pombe.momics #
momics ingest chroms -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/chrom_sizes_renamed/972hminus.chrom.sizes -g merged_yeasts S_pombe.momics
momics ingest seq -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/ref_genomes_renamed/972hminus.fa S_pombe.momics
momics ingest tracks -f H3=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/H3_concatenated_norm0999.bw S_pombe.momics
momics ingest tracks -f MNASE_15=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/MNASE_15min_concatenated_norm0999.bw S_pombe.momics

momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/C_glabrata/MNASE_25min_concatenated_norm0999.bw c_glabrata.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/C_glabrata/MNASE_35min_concatenated_norm0999.bw c_glabrata.momics

momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/MNASE_25min_concatenated_norm0999.bw S_pombe.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/MNASE_35min_concatenated_norm0999.bw S_pombe.momics


#-------------------------------

momics create S_eubayanus.momics #
momics ingest chroms -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/chrom_sizes_renamed/CDS216.chrom.sizes -g merged_yeasts S_eubayanus.momics
momics ingest seq -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/ref_genomes_renamed/CDS216.fa S_eubayanus.momics

momics ingest tracks -f H3=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_eubayanus/H3_concatenated_norm0999.bw S_eubayanus.momics
momics ingest tracks -f MNASE_15=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_eubayanus/MNASE_15min_concatenated_norm0999.bw S_eubayanus.momics
momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_eubayanus/MNASE_25min_concatenated_norm0999.bw S_eubayanus.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_eubayanus/MNASE_35min_concatenated_norm0999.bw S_eubayanus.momics

/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/momics_objects/Manon_momics/S_paradoxus.momics


momics create S_paradoxus.momics #
momics ingest chroms -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/chrom_sizes_renamed/CBS432.chrom.sizes -g merged_yeasts S_paradoxus.momics
momics ingest seq -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/ref_genomes_renamed/CBS432.fa S_paradoxus.momics
momics ingest tracks -f H3=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_paradoxus/H3_concatenated_norm0999.bw S_paradoxus.momics
momics ingest tracks -f MNASE_15=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_paradoxus/MNASE_15min_concatenated_norm0999.bw S_paradoxus.momics
momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_paradoxus/MNASE_25min_concatenated_norm0999.bw S_paradoxus.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_paradoxus/MNASE_35min_concatenated_norm0999.bw S_paradoxus.momics

momics create S_pombe.momics #
momics ingest chroms -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/chrom_sizes_renamed/972hminus.chrom.sizes -g merged_yeasts S_pombe.momics
momics ingest seq -f /pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/reference_genomes/5_yeasts/ref_genomes_renamed/972hminus.fa S_pombe.momics
momics ingest tracks -f H3=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/H3_concatenated_norm0999.bw S_pombe.momics
momics ingest tracks -f MNASE_15=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/MNASE_15min_concatenated_norm0999.bw S_pombe.momics
momics ingest tracks -f MNASE_25=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/MNASE_25min_concatenated_norm0999.bw S_pombe.momics
momics ingest tracks -f MNASE_35=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/Manon_MNase/normalized_bw_split_by_species/S_pombe/MNASE_35min_concatenated_norm0999.bw S_pombe.momics

/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/tracks_polymerase_normalized/RPA14.bw

momics create polymerase.momics #
momics ingest chroms -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/S288c.chrom.sizes -g S288c polymerase.momics
momics ingest seq -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/S288c.fa polymerase.momics
momics ingest tracks -f RPA14=/pasteur/helix/scratch/adecugis/all_momics/momics/inputs/tracks_polymerase_normalized/RPA14.bw polymerase.momics


cd /pasteur/appa/scratch/adecugis/all_momics/momics/momics_objects
momics delete -y glabrata.momics

momics create cerevisiae.momics #
momics ingest chroms -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/S288c.chrom.sizes -g S288c cerevisiae.momics
momics ingest seq -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/S288c.fa cerevisiae.momics

momics ingest tracks -f H3K4ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Cerevisae/H3K4ME3.bw cerevisiae.momics
momics ingest tracks -f H3K36ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Cerevisae/H3K36ME3.bw cerevisiae.momics
momics ingest tracks -f H3K27AC=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Cerevisae/H3K27AC.bw cerevisiae.momics
momics ingest tracks -f SCC1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Cerevisae/SCC1.bw cerevisiae.momics
momics ingest tracks -f RPB1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Cerevisae/STP1.bw cerevisiae.momics
momics ingest tracks -f SUA7=/pasteur/appa/scratch/adecugis/all_momics/momics/output/a_extracting_bw/cerevisiae/SUA7.bw cerevisiae.momics

"ASH1", "BRE1", "BRE2", "CBC2", "CCL1", "ABF1", "H3K4ME3_rescaled", "H3K27AC_rescaled", "H3K36ME3_rescaled", "MNase_rescaled", "RPB1_rescaled", "SCC1_rescaled", "STP1_rescaled", "CET1", "CSN12", "CTR9", "CYC8", "EAF7", "ELF1", "GLN3", "HDA1", "HMS2", "HOM6", "HPA3", "REB1", "SIR3", "SUA7"

momics ingest tracks -f MNASE=/pasteur/appa/scratch/adecugis/all_momics/momics/output/a_extracting_bw/cerevisiae/cerevisiae_MNASE_rescaled.bw cerevisiae.momics
momics ingest tracks -f STP1=/pasteur/appa/scratch/adecugis/all_momics/momics/output/a_extracting_bw/cerevisiae/cerevisiae_STP1_rescaled.bw cerevisiae.momics

momics ingest tracks -f ASH1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/ASH1.bw cerevisiae.momics
momics ingest tracks -f BRE1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/BRE1.bw cerevisiae.momics
momics ingest tracks -f BRE2=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/BRE2.bw cerevisiae.momics
momics ingest tracks -f CBC2=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/CBC2.bw cerevisiae.momics
momics ingest tracks -f CCL1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/CCL1.bw cerevisiae.momics
momics ingest tracks -f ABF1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/ABF1.bw cerevisiae.momics
momics ingest tracks -f CET1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/CET1.bw cerevisiae.momics
momics ingest tracks -f CSN12=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/CSN12.bw cerevisiae.momics
momics ingest tracks -f CTR9=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/CTR9.bw cerevisiae.momics
momics ingest tracks -f CYC8=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/CYC8.bw cerevisiae.momics
momics ingest tracks -f EAF7=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/EAF7.bw cerevisiae.momics
momics ingest tracks -f ELF1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/ELF1.bw cerevisiae.momics
momics ingest tracks -f GLN3=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/GLN3.bw cerevisiae.momics
momics ingest tracks -f HDA1=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/HDA1.bw cerevisiae.momics
momics ingest tracks -f HMS2=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/HMS2.bw cerevisiae.momics
momics ingest tracks -f HOM6=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/HOM6.bw cerevisiae.momics
momics ingest tracks -f HPA3=/pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracks_bw/HPA3.bw cerevisiae.momics

#---------------------

momics create glabrata.momics #
momics ingest chroms -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138.chrom.sizes -g CBS138 glabrata.momics
momics ingest seq -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138.fa glabrata.momics

momics ingest tracks -f H3K4ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/H3K4ME3.bw glabrata.momics
momics ingest tracks -f H3K36ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/H3K36ME3.bw  glabrata.momics
momics ingest tracks -f H3K27AC=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/H3K27AC.bw glabrata.momics
momics ingest tracks -f SCC1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/SCC1.bw glabrata.momics
momics ingest tracks -f RPB1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/RPB1.bw glabrata.momics
momics ingest tracks -f MNASE=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/MNASE.bw  glabrata.momics


momics ingest features -f blacklist=/pasteur/appa/scratch/adecugis/all_momics/momics/momics_objects/exclude.bed cerevisiae.momics 


momics ingest tracks -f H3K4ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Glabrata/H3K4_glabrata_SRR26031408_over_SRR26031394_CBS138_ratio.bw glabrata.momics
momics ingest tracks -f H3K27AC=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Glabrata/H3K27AC_Glabrata_SRR26031561_over_SRR26031448_CBS138_ratio.bw glabrata.momics
momics ingest tracks -f H3K36ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Glabrata/H3K36_Glab...ratio.bw glabrata.momics
momics ingest tracks -f SCC1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Glabrata/scc1_glabrata_LM184_over_LM188_CBS138_ratio.bw glabrata.momics
momics ingest tracks -f RPB1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Glabrata/Glabrata_RPB1_SRR26031000_over_SRR26030964_CBS138_ratio.bw glabrata.momics
momics ingest tracks -f MNASE=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Glabrata/Glabrata_Mnase_SRR24895593_resized.bw glabrata.momics
#---------------------------------
momics create glabrata.momics #
momics ingest chroms -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138.chrom.sizes -g CBS138 glabrata.momics
momics ingest seq -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/CBS138.fa glabrata.momics

momics ingest tracks -f H3K4ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/H3K4ME3.bw glabrata.momics
momics ingest tracks -f H3K36ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/H3K36ME3.bw glabrata.momics
momics ingest tracks -f H3K27AC=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/H3K27AC.bw  glabrata.momics
momics ingest tracks -f SCC1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/SCC1.bw glabrata.momics
momics ingest tracks -f MNASE=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/MNASE.bw  glabrata.momics
momics ingest tracks -f RPB1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Glabrata/RPB1.bw glabrata.momics

#---------------------------------
momics delete -y pombe.momics
momics create pombe.momics #
momics ingest chroms -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe.chrom.sizes -g Pombe pombe.momics
momics ingest seq -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe.fa pombe.momics

momics ingest tracks -f H3K4ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Pombe/GSM1201988_H3K4me3_ENA.bw pombe.momics
momics ingest tracks -f H3K36ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Pombe/GSM1201986_H3K36me3_ENA.bw pombe.momics
momics ingest tracks -f MNASE=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Pombe/Spombe_Mnase_ENA.bw  pombe.momics
momics ingest tracks -f RPB1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_built/Pombe/GSM6409929_01-13_Rpb1_ChIP-seq_WT_ENA_ENA.bw pombe.momics

momics ingest tracks -f H3K4ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Pombe/H3K4ME3.bw pombe.momics
momics ingest tracks -f H3K36ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Pombe/H3K36ME3.bw  pombe.momics
momics ingest tracks -f SCC1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Pombe/SCC1.bw pombe.momics
momics ingest tracks -f RPB1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Pombe/RPB1.bw  pombe.momics 
momics ingest tracks -f MNASE=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Pombe/MNASE.bw  pombe.momics

momics ingest tracks -f ABF1=/pasteur/appa/scratch/adecugis/handle_files/inputs/clean_input_dataset_scaled_90/Pombe/MNASE.bw  pombe.momics


momics ingest tracks -f H3K9ME3=/pasteur/appa/scratch/adecugis/handle_files/inputs/Pombe_H3k9Me3_ENA_norm.bw pombe.momics 


#making bw scc1 pombe
srun bowtie2 \
  -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  -U SRR24512031^unmapped_Pombe^T29NUU.1.gz \
  -S pombe.sam

srun samtools view -bS pombe.sam | samtools sort -o pombe.sorted.bam
srun samtools index pombe.sorted.bam

srun bamCoverage \
  -b pombe.sorted.bam \
  -o pombe.bw \
  --binSize 10 \
  --normalizeUsing CPM

#--------------------------------------------
#INGEST ALL TRACKS
momics create multi_tracks_cerevisiae.momics

momics ingest chroms \
  -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Cerevisiae.chrom.sizes \
  -g Cerevisiae \
  multi_tracks_cerevisiae.momics

momics ingest chroms -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/S288c.chrom.sizes -g S288c multi_tracks_cerevisiae.momics
momics ingest seq -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/S288c.fa multi_tracks_cerevisiae.momics


momics ingest seq \
  -f /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Cerevisiae.fa \
  multi_tracks_cerevisiae.momics

for bw in /pasteur/appa/scratch/adecugis/all_momics/momics/inputs/tracs_bw/*.bw; do
  name=$(basename "$bw" .bw)
  momics ingest tracks -f ${name}=${bw} multi_tracks_cerevisiae.momics
done

#---------------------
#MAKE MOMICS FOR EACH TRACK

srun --cpus-per-task=4 --mem=16G \
bowtie2 \
  -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  -U /pasteur/appa/scratch/adecugis/handle_files/outputs/output_SRR24512030/fastq/spikein/SRR24512030/SRR24512030^unmapped_Pombe^T29NUU.1.gz \
  -S SRR24512031_pombe.sam

srun --cpus-per-task=4 --mem=16G \
bowtie2 \
  -x /pasteur/appa/scratch/adecugis/handle_files/inputs/reference/Pombe \
  -U /pasteur/appa/scratch/adecugis/handle_files/outputs/output_SRR24512030/fastq/spikein/SRR24512031/SRR24512031^unmapped_Pombe^T29NUU.1.gz \
  -S SRR24512030_pombe.sam


-----------

srun --cpus-per-task=8 --mem=16G --time=02:00:00 bash -c "

samtools view -@ 8 -bS \
/pasteur/appa/scratch/adecugis/handle_files/inputs/scratch_inputs/SRR24512031_pombe.sam \
| samtools sort -@ 8 -o SRR24512031_pombe.sorted.bam

# Index
samtools index SRR24512031_pombe.sorted.bam

# IP / Input ratio (no log)
bamCompare \
  -b1 SRR24512030_pombe.sorted.bam \
  -b2 SRR24512031_pombe.sorted.bam \
  --operation ratio \
  --binSize 10 \
  --normalizeUsing CPM \
  --pseudocount 1 \
  -p 8 \
  -o SRR24512030_over_31_ratio.bw
"

