#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 00-20:30:00
#SBATCH --mem=600G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J Part1_Tele04
#SBATCH -o Part1_Tele04_%j.out
#SBATCH -e Part1_Tele04_%j.err

#############################################
# Metabarcoding pipeline – OBITools 1.2.11 + Sumaclust
# Author: turvois
# Date: $(date)
#############################################

#=============
# Variables
#=============
WORKDIR="/work/user/turvois/metabar_Tele04"
OUTDIR="$WORKDIR/output"

BarcodeName="12STeleo4"
Fwd="GTGCCAGCCACCGCGGTT"
Rev="GTGGGGTATCTAATCCCAGTTTG"
RevC="CAAACTGGGATTAGATACCCCAC"
MinSeqLength=110
MinCount=100
ecoPCRdb="/work/user/equemere/metabarcoding/ecoPCR/ecoPCR_12sTELEO4.fasta"
input_fastq="$WORKDIR/merged_all.fastq"

NCPU=$SLURM_CPUS_PER_TASK

#=============
# Modules
#=============
module load bioinfo/Cutadapt/4.3
module load devel/Miniconda/Miniconda3
module load bioinfo/OBITools/1.2.11
module load bioinfo/Sumaclust/1.0.36

#=============
# Logs
#=============
mkdir -p "$OUTDIR/logs" "$OUTDIR/stats"

#=============
# Pipeline
#=============

echo ">>> Step 1: Primer trimming with Cutadapt"
cutadapt -g $Fwd -O 10 -e 0.1 --discard-untrimmed -j $NCPU \
  $input_fastq -o $OUTDIR/cutadapt_fwd.fastq --revcomp \
  2> $OUTDIR/logs/cutadapt_fwd.log

cutadapt -a $RevC -O 10 -e 0.1 --discard-untrimmed -j $NCPU \
  $OUTDIR/cutadapt_fwd.fastq -o $OUTDIR/${BarcodeName}_trim.fastq \
  --revcomp --match-read-wildcards \
  2> $OUTDIR/logs/cutadapt_rev.log

obistat $OUTDIR/${BarcodeName}_trim.fastq > $OUTDIR/stats/step1_cutadapt.stats


echo ">>> Step 2: Dereplication with obiuniq (obitools)"
obiuniq -m sample $OUTDIR/${BarcodeName}_trim.fastq \
  > $OUTDIR/${BarcodeName}_uniq.fasta

obiannotate -k count -k merged_sample $OUTDIR/${BarcodeName}_uniq.fasta \
  > $OUTDIR/${BarcodeName}_uniq_annot.fasta

obistat $OUTDIR/${BarcodeName}_uniq_annot.fasta > $OUTDIR/stats/step2_derep.stats


echo ">>> Step 3: Filtering (min length $MinSeqLength, min count $MinCount)"
obigrep -l $MinSeqLength -p "count>=$MinCount" \
  $OUTDIR/${BarcodeName}_uniq_annot.fasta \
  > $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}.fasta

obistat $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}.fasta \
  > $OUTDIR/stats/step3_filter.stats


#------------------------------------------------------
echo ">>> Step 4: Clustering with Sumaclust at 97%"
#------------------------------------------------------
sumaclust \
  $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}.fasta \
  -t 0.97 \
  > $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97.fasta \
  2> $OUTDIR/logs/sumaclust_c${MinCount}_l${MinSeqLength}.log

obiselect -c cluster -n 1 --merge sample -M -f count \
  $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97.fasta \
  > $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97_agg.fasta

obistat $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97_agg.fasta > $OUTDIR/stats/step4_c${MinCount}_l${MinSeqLength}_cluster97.stats


#------------------------------------------------------
echo ">>> Step 5: Sorting & OTU tables"
#------------------------------------------------------
obisort -k count -r \
  $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97_agg.fasta \
  > $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97_agg_sort.fasta

obitab -o $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97_agg_sort.fasta \
  > $OUTDIR/${BarcodeName}_c${MinCount}_l${MinSeqLength}_cluster97_sort_agg.tab


#=============
# Done
#=============
echo ">>> Pipeline finished successfully!"