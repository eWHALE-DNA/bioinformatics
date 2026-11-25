#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 00-20:30:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J Part2_Tele04
#SBATCH -o Part2_Tele04_%j.out
#SBATCH -e Part2_Tele04_%j.err

#############################################
# Metabarcoding pipeline – OTU FASTA splitter
# Author: turvois
# Date: $(date)
#############################################

#=============
# Variables
#=============

WORKDIR="/work/user/turvois/metabar_Tele04"
OUTDIR="$WORKDIR/output"
INPUT_FASTA="$OUTDIR/OTU_1000_merged_1perc.fasta"
SINGLE_OTU_DIR="$OUTDIR/single_OTU"

#=============
# Setup
#=============

mkdir -p $OUTDIR/single_OTU

#=============
# Pipeline
#=============

echo ">>> Splitting OTUs into individual FASTA files"

awk '
/^>/ {
    # Keep only the part before first space or semicolon
    split($0, a, /[ ;]/)
    name = a[1]
    sub(/^>/, "", name)                 # remove leading >
    gsub(/[^a-zA-Z0-9]/, "_", name)     # sanitize filename
    file = "'"$SINGLE_OTU_DIR"'/" name ".fasta"
}
{
    print > file
}
' "$INPUT_FASTA"

echo ">>> OTU splitting completed successfully"