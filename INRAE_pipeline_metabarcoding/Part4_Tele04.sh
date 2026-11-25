#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 00-20:30:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J Part4_Tele04
#SBATCH -o Part4_Tele04_%j.out
#SBATCH -e Part4_Tele04_%j.err

#############################################
# Metabarcoding pipeline – Concatenate assigned OTUs + generate OTU tables
# Author: turvois
# Date: $(date)
#############################################

#=============
# Variables
#=============
WORKDIR="/work/user/turvois/metabar_Tele04"
OUTDIR="$WORKDIR/output"
SINGLE_OTU_DIR="$OUTDIR/single_OTU"

#=============
# Modules
#=============
module load bioinfo/Cutadapt/4.3
module load devel/Miniconda/Miniconda3
module load bioinfo/OBITools/1.2.11
module load bioinfo/Sumaclust/1.0.36

#------------------------------------------------------
echo ">>> Step 1: Concatenate all assigned OTUs + any failed OTU"
#------------------------------------------------------
cat "$SINGLE_OTU_DIR"/*_assigned.fasta \
    "$SINGLE_OTU_DIR/VH00807_358_AAGJFLJM5_1_1101_23155_1189.fasta" \
    > "$OUTDIR/OTU_1000_merged_1perc_assigned.fasta"

#------------------------------------------------------
echo ">>> Step 2: Sort sequences & generate OTU table"
#------------------------------------------------------
# Sort FASTA by count
obisort -k count -r \
  "$OUTDIR/OTU_1000_merged_1perc_assigned.fasta" \
  > "$OUTDIR/OTU_1000_merged_1perc_assigned_sort.fasta"

# Generate OTU table (counts per OTU)
obitab -o "$OUTDIR/OTU_1000_merged_1perc_assigned_sort.fasta" \
  > "$OUTDIR/OTU_1000_merged_1perc_assigned_sort.tab"

#=============
# Done
#=============
echo ">>> Pipeline finished successfully!"