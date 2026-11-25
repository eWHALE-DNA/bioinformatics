#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 00-20:30:00
#SBATCH --mem=600G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J Part0_Tele04
#SBATCH -o Part0_Tele04_%j.out
#SBATCH -e Part0_Tele04_%j.err


#############################################
# Metabarcoding pipeline – FASTQ Paired-End Merging with VSEARCH
# Author: turvois
# Date: $(date)
#############################################

#=============
# Variables
#=============
WORKDIR="/work/user/turvois/metabar_Tele04"
DATADIR="$WORKDIR/raw_fastq"
OUTDIR="$WORKDIR/merged_data"


#=============
# Modules
#=============
module load bioinfo/VSEARCH/2.22.1

#=============
# Pipeline
#=============

# Create output directories
mkdir -p "$OUTDIR"

# Move to data directory
cd "$DATADIR" || exit 1

# Count total number of R1 files
R1_count=$(ls -1 *_R1.fastq 2>/dev/null | wc -l)
echo "Total number of R1 files: $R1_count"

# Initialize counters
merged_count=0
total_count=0
failed_count=0

# Loop over all R1 files
for R1_file in *_R1.fastq; do
    R2_file=${R1_file/_R1.fastq/_R2.fastq}
    total_count=$((total_count + 1))

    if [[ -f "$R2_file" ]]; then
        full_base_name=$(basename "$R1_file" _R1.fastq)
        sample_name=$(echo "$full_base_name" | cut -d'_' -f1)

        vsearch --fastq_mergepairs "$R1_file" \
                --reverse "$R2_file" \
                --fastqout "$OUTDIR/${full_base_name}_merged.fastq" \
                --fastq_minovlen 10 \
                --fastq_maxdiffs 10 \
                --fastq_allowmergestagger \
                --fastq_qmax 50 \
                --label_suffix ";sample=${sample_name};" \
                --quiet

        if [[ $? -eq 0 ]]; then
            merged_count=$((merged_count + 1))
        else
            failed_count=$((failed_count + 1))
            echo "??  Merge failed for sample: $full_base_name"
        fi
    else
        echo "??  Missing R2 file for $R1_file"
        failed_count=$((failed_count + 1))
    fi

    # Progress report every 50 files
    if [[ $((total_count % 50)) -eq 0 ]]; then
        echo "? Progress: $total_count processed | $merged_count merged | $failed_count failed"
    fi
done

# Concatenate all merged FASTQ files
cat "$OUTDIR"/*_merged.fastq > "$WORKDIR/merged_all.fastq"

if [[ $? -eq 0 ]]; then
    echo "? merged_all.fastq successfully created"

    seq_count=$(wc -l < "$WORKDIR/merged_all.fastq")
    seq_count=$((seq_count / 4))
    echo "?? Total merged sequences: $seq_count"

    file_size=$(du -h "$WORKDIR/merged_all.fastq" | cut -f1)
    echo "?? Output file size: $file_size"
else
    echo "? ERROR: Failed to concatenate merged FASTQ files"
    exit 1
fi