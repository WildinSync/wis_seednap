## Bash script to demultiplex, put in correct order, filter out sequences which have nucleotides between the tag and the primer,
# but! not cutting the primers (it will be done later to be able to group all samples - demultiplexed and non demultiplexed)

# Use:
# bash scripts/demultiplex_ligation.sh [path_to_folder_with_libraries] [path_to_all_samples] [path_to_config.sh]
# bash scripts/demultiplex_ligation.sh [path_to_folder_with_libraries] [path_to_metadata_lab]] [path_to_config.sh]
# Metadata lab is more efficient than all_samples which is annoying to work with
# Metadata lab is formatted along ELE's group standards - it needs at least the columns eventID (sample_name) and library (full name of the pooled library)
# It creates a folder raw_for_dada2 in the folder of raw_path_ligation [first argument] or just locally for now in raw/ligation

# Examples
# bash scripts/demultiplex_ligation.sh /home/shared/edna/ma_fr_spygen_2020/  /home/shared/edna/ma_fr_spygen_2020/all_samples.csv config_teleo.sh
# bash scripts/demultiplex_ligation.sh /home/shared/edna/raw/MA_Calanques_2023/ /home/shared/edna/raw/MA_Calanques_2023/metadata/metadata_lab_ma_fr_calanq_2023.csv config_teleo.sh

# For debug
# raw_path_ligation="/home/shared/edna/raw/MA_Calanques_2023/"
# metadata_lab="/home/shared/edna/raw/MA_Calanques_2023/metadata/metadata_lab_ma_fr_calanq_2023.csv"
# source config/config_teleo.sh 

# Load variables
raw_path_ligation=$1
metadata_lab=$2
source $3

mkdir -p ${raw_path_ligation}/raw_for_dada2/

# For primer trimming 
primer_fwd=${primer_F}
primer_fwd_rc=$(echo "${primer_F}" | tr 'ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv' 'TACGAtacgaNnRrYySsWwMmKkVvHhDdBb' | rev)

primer_rev=${primer_R}
primer_rev_rc=$(echo "${primer_R}" | tr 'ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv' 'TACGAtacgaNnRrYySsWwMmKkVvHhDdBb' | rev)

# Execute to obtain the fasta demultiplexing files to demultiplex ligation-based libraries
Rscript scripts/generate_cutadapt_ligation.R $metadata_lab

# ---------------------------------------------------------------------- # 
# DEMULTIPLEX
# ---------------------------------------------------------------------- # 

path_cutadapt_files="outputs/00_demultiplex_ligation/cutadapt_tags/"
mkdir -p outputs/00_demultiplex_ligation/demultiplex
output_dir_tag_trim="outputs/00_demultiplex_ligation/demultiplex"

# Run command
for tag_file in `ls ${path_cutadapt_files}` 
do
  
  # Get run prefix
  run_prefix=`echo $tag_file | cut -d '.' -f 1` 
  # Get R1
  r1_in=$(ls "${raw_path}" | grep "^${run_prefix}.*_R1\.fastq\.gz$")
  # Get R2
  r2_in=$(ls "${raw_path}" | grep "^${run_prefix}.*_R2\.fastq\.gz$")

  # Execute cutadapt command
  cutadapt \
    -j $CORES \
    --discard-untrimmed \
    -e 0.0 --no-indels \
    -a file:"${path_cutadapt_files}/${tag_file}" -A file:"${path_cutadapt_files}/${tag_file}" \
    -o "$output_dir_tag_trim"/{name}.R1.fastq.gz -p "$output_dir_tag_trim"/{name}.R2.fastq.gz \
    "${raw_path}$r1_in" "${raw_path}$r2_in"
done

# ---------------------------------------------------------------------- # 
# PRIMER TRIMING
# ---------------------------------------------------------------------- # 

# Run primmer trimming 
# Here we want to detected the primers (cut where we dont have primers, but not trimming the primmers)
# We also want to make sure the primers are anchored, to avoid retaining reads with nucleotides between tag and primer
# Which in the past was source of issues 
# We also want to re-order the reads to make sure all R1 and R2 reads are in the expected order (ligation-induced)
# End file: R1/R2 in correct order, primers kept but only reads with anchored primers detected are retained and written in output directory

mkdir -p outputs/00_demultiplex_ligation/primer_detection
output_dir_primers_detect="outputs/00_demultiplex_ligation/primer_detection"

samples=$(ls "$output_dir_tag_trim"| grep -v "unknown" |cut -d. -f1 | sort | uniq)
echo $samples

# Proceed with the expected order of reads in R1/R2 (round 1 of trim)
for s in $samples;
do
 cutadapt \
	   -j $CORES \
       --action=none \
       --discard-untrimmed \
	   -e 0.1 --no-indels -m 20 \
	   -g "^${primer_fwd}...${primer_rev_rc}" -A "^${primer_rev}...${primer_fwd_rc}" \
	   -o $output_dir_primers_detect/trim_round1_${s}.R1.fastq.gz -p $output_dir_primers_detect/trim_round1_${s}.R2.fastq.gz \
	   $output_dir_tag_trim/${s}.R1.fastq.gz $output_dir_tag_trim/${s}.R2.fastq.gz ;
done

# Proceed with the unexpected order of reads in R1/R2 (round 2 of trim)
for s in $samples;
do
 cutadapt \
	   -j $CORES \
       --action=none \
       --discard-untrimmed \
	   -e 0.1 --no-indels -m 20 \
	   -A "^${primer_fwd}...${primer_rev_rc}" -g "^${primer_rev}...${primer_fwd_rc}"\
	   -o $output_dir_primers_detect/trim_round2_${s}.R1.fastq.gz -p $output_dir_primers_detect/trim_round2_${s}.R2.fastq.gz \
     $output_dir_tag_trim/${s}.R1.fastq.gz $output_dir_tag_trim/${s}.R2.fastq.gz ;
done

# ---------------------------------------------------------------------- # 
# MERGING
# ---------------------------------------------------------------------- # 

# Now we merge the R1 and R2 back together in the correct orientation 
# Note that half of the reads now on R1 are originally from R2 and vice-versa 
mkdir -p outputs/00_demultiplex_ligation/realigned
output_dir_realigned="outputs/00_demultiplex_ligation/realigned"

samples_c=$(ls "$output_dir_primers_detect"| grep -v "unknown" | cut -d_ -f3-| cut -d. -f1 | sort | uniq)
echo $samples_c

for s in $samples_c;
do
	cat "$output_dir_primers_detect"/trim_round1_${s}.R1.fastq.gz "$output_dir_primers_detect"/trim_round2_${s}.R2.fastq.gz > $output_dir_realigned/${s}.R1.fastq.gz
	cat "$output_dir_primers_detect"/trim_round1_${s}.R2.fastq.gz "$output_dir_primers_detect"/trim_round2_${s}.R1.fastq.gz > $output_dir_realigned/${s}.R2.fastq.gz
done

# Unzipp all merged
gunzip ${output_dir_realigned}/*.gz

# Now copy it back to raw path in a sub-directory
cp ${output_dir_realigned}/*.fastq ${raw_path_ligation}/raw_for_dada2/