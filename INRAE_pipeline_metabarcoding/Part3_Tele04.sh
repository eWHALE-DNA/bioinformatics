#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 00-20:30:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J Part3_Tele04
#SBATCH -o Part3_Tele04_%j.out
#SBATCH -e Part3_Tele04_%j.err

#############################################
# Metabarcoding pipeline – Taxonomic assignment (ecotag)
# Author: turvois
# Date: $(date)
#############################################

#=============
# Variables
#=============
WORKDIR="/work/user/turvois/metabar_Tele04"
OUTDIR="$WORKDIR/output"
SINGLE_OTU_DIR="$OUTDIR/single_OTU"

ECO_DATABASE="/work/bank2/ecopcr/EMBL-ENA-20211221-ECOPCR/embl_20211221"
ecoPCRdb="/work/user/equemere/metabarcoding/ecoPCR/ecoPCR_12sTELEO4.fasta" 
  # FASTA database generated with ecoPCR using TELEO4 primers

NCPU=8

#=============
# Modules
#=============
module load bioinfo/Cutadapt/4.3
module load devel/Miniconda/Miniconda3
module load bioinfo/OBITools/1.2.11
module load bioinfo/Sumaclust/1.0.36

#=============
# Setup
#=============
mkdir -p "$OUTDIR/logs"

#=============
# Pipeline
#=============
echo ">>> Starting taxonomic assignment with ecotag"

ls "$SINGLE_OTU_DIR"/*.fasta | parallel -j "$NCPU" \
  ecotag \
    -d "$ECO_DATABASE" \
    -m 0.95 \
    -R "$ecoPCRdb" \
    {} \
    '>' {.}_assigned.fasta \
    2> "$OUTDIR/logs/{.}_assigned.log"

echo ">>> ecotag assignment completed successfully"