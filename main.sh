# Execute primer trimming for one marker 
# V. Marques
# Last updated: 07/07/2025
# Usage: bash trimming.sh config.sh

# Source functions 
source scripts/bash_functions.sh

# Config source
# source config/config_amph.sh
source $1

# Create folder sctructure 
mkdir -p outputs/01_trim/$marker
mkdir -p outputs/02_dada2/$marker
mkdir -p outputs/03_taxo/$marker
output_trim="outputs/01_trim/"${marker}

# For primer trimming 
primer_fwd=$primer_F
primer_fwd_rc=$(reverse_complement "$primer_fwd")

primer_rev=$primer_R
primer_rev_rc=$(reverse_complement "$primer_rev")

# Remove the I files 
find $raw_path -type f -name '*_I[12]_001.fastq.gz' | xargs rm -v

# ---------------------------------------------------------------------- # 
# DEMULTIPLEXING IF NECESSARY
# ---------------------------------------------------------------------- # 

# TBA

# ---------------------------------------------------------------------- # 
# PRIMER TRIMING
# ---------------------------------------------------------------------- # 

# Run primmer trimming 
#### now primer trimming, only if we dont do it here and not in DADA2

samples=$(ls "$raw_path"| grep -v "unknown" | grep -v "Undetermined" |cut -d_ -f1 | sort | uniq)
echo $samples

# For loop is not working properly
# So here is to while loop

# Normal direction
while read -r s; do
  # Print
  echo $s
  # Get raw reads
  r1_in=$(ls "$raw_path" | grep "^${s}_R1.fastq.gz$")
  r2_in=$(ls "$raw_path" | grep "^${s}_R2.fastq.gz$")

  # Pass 1: Trim with -g / -g 
  # Run command
    cutadapt \
    -j $CORES \
    -e 0.1 \
    -m 20 \
    -g "${primer_fwd}" -G "${primer_rev}" \
    --untrimmed-output $output_trim/untrimmed_${s}.R1.fastq \
    --untrimmed-paired-output $output_trim/untrimmed_${s}.R2.fastq \
    -o $output_trim/${s}.R1_TEMPORARY.fastq -p $output_trim/${s}.R2_TEMPORARY.fastq \
    $raw_path/$r1_in $raw_path/$r2_in >> "logs/${marker}_log_trim_${s}.txt"

  # Pass 2: Trim with -a / -A - do not discard untrimmed
  cutadapt \
    -j $CORES \
    -e 0.1 \
    -m 20 \
    -a "${primer_rev_rc}" -A "${primer_fwd_rc}" \
    -o $output_trim/${s}.R1.fastq -p $output_trim/${s}.R2.fastq \
    $output_trim/${s}.R1_TEMPORARY.fastq $output_trim/${s}.R2_TEMPORARY.fastq >> "logs/${marker}_log_trim_${s}_2nd_trim.txt"

done <<< "$samples"  

# For now, remove the untrimmed (to investigate later and save elsewhere in complete version)
rm outputs/01_trim/$marker/untrimmed*
rm outputs/01_trim/$marker/*TEMPORARY*

# ------ # 
# And then do DADA2 on them to finish the processing

# Apply dada2 and taxo assignment on sequences
# Right now, dada2 is run all together, in further version run it by batch of sequencers and pool them later
Rscript scripts/dada2_process.R $marker
echo ">>${method_assignment}<<"
# Now apply taxonomic assignment depending on the method chosen 

if [[ "$method_assignment" == "ecotag" ]]; then
    echo "Execute ecotag assignment method"
    # Activate conda for obitools (not very clean but python conflict due to old obitools)
    conda activate ref_database
    # Run scripts
    bash scripts/taxo_ecotag_marker.sh "$path_ecotag_tree" "$path_ecotag_fasta" "$marker"
    # Link assignment with abundance table
    Rscript scripts/link_taxo_files_ecotag.R "$marker"
    # Re-activate the original conda
    conda activate metabarcoding

elif [[ "$method_assignment" == "dada2" ]]; then
    # Check if the current user is in the sudo group
    echo "Execute dada2 assignment method"
    Rscript scripts/taxo_dada2_marker.R "$marker" "$path_dada_all" "$path_dada_species"

elif [[ "$method_assignment" == "decipher" ]]; then
    echo "Execute decipher assignment method"
    Rscript scripts/taxo_decipher_marker.R "$marker" "$path_decipher_trained"

elif [[ "$method_assignment" == "blast" ]]; then
    echo "Execute blast assignment method"
    # Run blastn script and table creation
    bash scripts/taxo_blast_marker.sh "$path_blast_fasta" "$marker"

else
    echo "Unknown assignment method: $method_assignment"
fi

# Now the formatting - only for dada2 for now
if [[ ${method_assignment} == "dada2" && ${convert_for_GBIF} == "TRUE" ]]; then 
    Rscript scripts/RDP_to_GBIFinput.R outputs/${marker}_dada2RDP.csv
elif [[ ${method_assignment} == "ecotag" && ${convert_for_GBIF} == "TRUE" ]]; then 
    Rscript scripts/ecotag_to_GBIFinput.R outputs/${marker}_ecotag.csv
fi

# DONE ! 
