#!/bin/sh

# The following implements the OBITools v1.2 metabarcoding pipeline.
# General idea is that users do not need to edit this script.
# This script takes user-defined parameters from the "config.ini" file and will output a .csv file containing clustered reads with taxon assignments.
# Please see "README.txt" for implementations details.


# Load configuration from config.ini
CONFIG_FILE="config.ini"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found!"
    exit 1
fi

# Read values from config.ini
DataDIR=$(awk -F "=" '/^DataDIR/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
InputFASTA=$(awk -F "=" '/^InputFASTA/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
Label=$(awk -F "=" '/^Label/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
Threads=$(awk -F "=" '/^Threads/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
SwarmD=$(awk -F "=" '/^SwarmD/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
Marker=$(awk -F "=" '/^Marker/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
# Read Ecotag database paths from config.ini
RefDB_12S=$(awk -F "=" '/^RefDB_12S/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
RefDB_CO1=$(awk -F "=" '/^RefDB_CO1/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
RefFASTA_12S=$(awk -F "=" '/^RefFASTA_12S/ {print $2}' "$CONFIG_FILE" | tr -d ' ')
RefFASTA_CO1=$(awk -F "=" '/^RefFASTA_CO1/ {print $2}' "$CONFIG_FILE" | tr -d ' ')

# Check if required parameters are set
if [[ -z "$DataDIR" || -z "$InputFASTA" || -z "$Label" || -z "$Threads" || -z "$SwarmD" || -z "$Marker" || -z "$RefDB_12S" || -z "$RefDB_CO1" || -z "$RefFASTA_12S" || -z "$RefFASTA_CO1" ]]; then
    echo "Error: One or more required parameters are missing from config.ini"
    exit 1
fi

echo "Using Data Directory: $DataDIR"
echo "Using Input FASTA: $InputFASTA"
echo "Using Label: $Label"
echo "Using Threads: $Threads"
echo "Using Swarm d-value: $SwarmD"
echo "Using Metabarcoding Marker: $Marker"
echo "Local Reference Database for 12S: $RefDB_12S"
echo "Local Reference Database for CO1: $RefDB_CO1"
echo "Fasta file for 12S: $RefFASTA_12S"
echo "Fasta file for CO1: $RefFASTA_CO1"

# Create res/ directory if it doesn't exist
mkdir -p "$DataDIR/res"

# Copy and rename the input FASTA file
cp "$DataDIR/$InputFASTA" "$DataDIR/res/${Label}.filtered.fasta"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to copy input FASTA file."
    exit 1
fi

# Change to working directory
cd "$DataDIR/res" || exit 1

#Load the paths for necessary packages and R scripts here

echo "Filtering input fasta file based on read length thresholds for marker ${Marker}..."

# Filter sequences based on the specified marker
if [[ "$Marker" == "12S" ]]; then
    obigrep -p 'seq_length>140' -p 'seq_length<190' -s '^[ACGT]+$' "$DataDIR/res/${Label}.filtered.fasta" > "$DataDIR/res/${Label}.filtered_length.fasta"
elif [[ "$Marker" == "CO1" ]]; then
    obigrep -p 'seq_length>290' -p 'seq_length<340' -s '^[ACGT]+$' "$DataDIR/res/${Label}.filtered.fasta" > "$DataDIR/res/${Label}.filtered_length.fasta"
else
    echo "Error: Invalid Marker specified. Choose either 12S or CO1."
    exit 1
fi

echo "Done"

echo "Creating plate statistics..."
#Create the plate statistics
obistat -c sample -a seq_length "$DataDIR/res/${Label}.filtered_length.fasta" > "$DataDIR/res/sample_stats.txt" &&
echo "Done"

echo "Collapsing identical sequences..."
#Collapse identical sequences, retain info on frequency of each unique sequence
obiuniq -m sample "$DataDIR/res/${Label}.filtered_length.fasta" > "$DataDIR/res/${Label}.unique.fasta" &&
echo "Done"

#Exchange the identifier to a short index
obiannotate --seq-rank "$DataDIR/res/${Label}.unique.fasta" | obiannotate --set-identifier '"'${Label}'_%09d" % seq_rank' > "$DataDIR/res/${Label}.new.fasta" &&
#obiannotate --seq-rank "$DataDIR/res/${Label}.unique.fasta" | obiannotate --set-identifier "${Label}_%09d" % seq_rank' > "$DataDIR/res/${Label}.new.fasta" &&

echo "Getting .tab file..."
#Get the tab file. This file contains sequence identifier and read abundance in every sample
obitab -o "$DataDIR/res/${Label}.new.fasta" > "$DataDIR/res/${Label}.new.tab" &&
echo "Done"

echo "Removing chimeras..."
#Remove chimeras
##Convert into vsearch format
owi_obifasta2vsearch_gg -i "$DataDIR/res/${Label}.new.fasta" -o "$DataDIR/res/${Label}.vsearch.fasta" &&

##Remove chimeras by runing UCHIME in Vsearch
vsearch --uchime_denovo "$DataDIR/res/${Label}.vsearch.fasta" --sizeout --nonchimeras "$DataDIR/res/${Label}.nonchimeras.fasta" --chimeras "$DataDIR/res/${Label}.chimeras.fasta" \
	--uchimeout "$DataDIR/res/${Label}.uchimeout.txt" &&
echo "Done"

echo "Converting to single line fasta for SWARM"
##SWARM functions best with single line fasta files
awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "$DataDIR/res/${Label}.nonchimeras.fasta" > "$DataDIR/res/${Label}.nonchimeras_singleline.fasta" &&
echo "Done"

echo "Clustering with SWARM and d=${SwarmD}"
#Clustering using SWARM with user-specified d and "fastidious" if d=1
if [[ "$SwarmD" == 1 ]]; then
        swarm -d "$SwarmD" -f -t "$Threads" -z -o "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_output" -s "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_stats" \
        -w "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_seeds.fasta" "$DataDIR/res/${Label}.nonchimeras_singleline.fasta"
else
    	swarm -d "$SwarmD" -t "$Threads" -z -o "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_output" -s "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_stats" \
        -w "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_seeds.fasta" "$DataDIR/res/${Label}.nonchimeras_singleline.fasta"
fi

#Recount after SWARM
owi_recount_swarm "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_output" "$DataDIR/res/${Label}.new.tab" &&

##Removing singletons
# Edit header ";size=" to "; size="
sed -i 's/;size=/; size=/g' "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_seeds.fasta" &&
obigrep -p 'size>1' "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_seeds.fasta" > "$DataDIR/res/${Label}_${SwarmD}nc_seeds_nonsingletons.fasta" &&

#LULU, identifying redundant MOTUs and pseudogenes
## Remove the ";size=N;" part of the identifier of the MOTU to get a cleaner file.
sed -i 's/ size=\([0-9]\+\);//g' "$DataDIR/res/${Label}_${SwarmD}nc_seeds_nonsingletons.fasta" &&

##More stringent clustering step to reduce MOTUs
vsearch --threads "$Threads" --usearch_global "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_seeds.fasta" --db "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_seeds.fasta" \
	--self --id .84 --iddef 1 --userout "$DataDIR/res/match_list_${SwarmD}nc.txt" -userfields query+target+id \
	--maxaccepts 0 --query_cov .9 --maxhits 10 &&

# Get first and last sample columns
first_sample_col=$(awk -F";" 'NR==1 { for (i=1; i<=NF; i++) if ($i ~ /^sample\./) {print i; exit} }' "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_output.counts.csv")
last_sample_col=$(awk -F";" 'NR==1 { for (i=NF; i>=1; i--) if ($i ~ /^sample\./) {print i; exit} }' "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_output.counts.csv")

if [[ -z "$first_sample_col" || -z "$last_sample_col" ]]; then
    echo "Error: Could not determine first and last sample columns in CSV file."
    exit 1
fi

# Process the counts file for LULU
cut -d ";" -f1,"${first_sample_col}"-"${last_sample_col}" "$DataDIR/res/${Label}_SWARM_${SwarmD}nc_output.counts.csv" > "$DataDIR/res/${Label}_${SwarmD}nc.counts_LULU.txt" &&
sed -i 's/;/\t/g' "$DataDIR/res/${Label}_${SwarmD}nc.counts_LULU.txt" &&

Rscript <<EOF
# check if pointing to R packages in correct location with: Rscript -e '.libPaths()'
# load necessary package
library(lulu)
# read input files
matchlist <- read.table("$DataDIR/res/match_list_${SwarmD}nc.txt", header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
otutab <- read.csv("$DataDIR/res/${Label}_${SwarmD}nc.counts_LULU.txt", sep='\t', header=TRUE, as.is=TRUE, row.names=1)
# run LULU curation
curated_result <- lulu(otutab, matchlist)
# save output
write.csv(curated_result$curated_table,"$DataDIR/res/${Label}_${SwarmD}nc.LULU.curated.csv")
EOF


# Set ecotag database paths based on marker type
if [[ "$Marker" == "12S" ]]; then
    RefDB="$RefDB_12S"
    RefFASTA="$RefFASTA_12S"
elif [[ "$Marker" == "CO1" ]]; then
    RefDB="$RefDB_CO1"
    RefFASTA="$RefFASTA_CO1"
else
    echo "Error: Invalid Marker specified for ecotag. Choose either 12S or CO1."
    exit 1
fi

echo "ready for ecotagging"
exit 1

# Run ecotag with the correct database
#ECOTAG_INPUT="$DataDIR/res/${Label}_${SwarmD}nc_seeds_nonsingletons.fasta"
#ECOTAG_OUTPUT="$DataDIR/res/${Label}_${SwarmD}nc_ecotag.fasta"

ecotag -d "$RefDB" -R "$RefFASTA" "$ECOTAG_INPUT" > "$ECOTAG_OUTPUT" &&

#Add High level taxa
owi_add_taxonomy "$DataDIR/res/${Label}_${SwarmD}nc_ecotag.fasta" &&

# Change the commas for dots in the column referring to the identity to the reference sequence
sed -i -e 's/,/./g' "$DataDIR/res/${Label}_${SwarmD}nc_ecotag.fasta.annotated.csv" &&

#Combine ecotag and abundance files
## Format LULU results
sed -i -e 's/""/"id"/g' "$DataDIR/res/${Label}_${SwarmD}nc.LULU.curated.csv" &&
sed -i -e 's/"//g' "$DataDIR/res/${Label}_${SwarmD}nc.LULU.curated.csv" &&
sed -i -e 's/,/;/g' "$DataDIR/res/${Label}_${SwarmD}nc.LULU.curated.csv" &&

## Combine files
owi_combine -i "$DataDIR/res/${Label}_${SwarmD}nc_ecotag.fasta.annotated.csv" -a "$DataDIR/res/${Label}_${SwarmD}nc.LULU.curated.csv" -o "$DataDIR/res/${Label}_${SwarmD}nc.combMOTUs.csv" &&

# In R. Edit the file adding the sequence from the taxonomy table
Rscript <<EOF
table_abundances <- read.csv("$DataDIR/res/${Label}_${SwarmD}nc.combMOTUs.csv", sep=";", header=T)
taxonomy_sequences <- read.csv("$DataDIR/res/${Label}_${SwarmD}nc_ecotag.fasta.annotated.csv", sep=";", header=T)
just_sequences <- taxonomy_sequences[,c("id", "sequence")]
table_abundance_sequences <- merge(table_abundances, just_sequences, by="id")
write.csv2(table_abundance_sequences, "$DataDIR/res/${Label}_${SwarmD}nc.combMOTUs.seqs.csv", row.names=F)
EOF

#Final step
refine_MOTU_table -i "$DataDIR/res/${Label}_${SwarmD}nc.combMOTUs.seqs.csv" -o "$DataDIR/res/${Label}_${SwarmD}nc.final.csv" -t 0.00 -r 5 -a yes -s 0.96 &&

echo "Pipeline completed successfully."
