# Config parameters
marker="teleo" # No majuscules here
raw_path="/home/shared/edna/raw/ma_dj_djibou_2023/"
primer_F="ACACCGCCCGTCACTCT"
primer_R="CTTCCGGTACACTTACCATG"
method_demultiplex="" # among: primer_trim / ligation_trim
convert_for_GBIF="TRUE" # values= TRUE or FALSE

# Assignment parameters
method_assignment="dada2" # among: "dada2 / decipher / ecotag / blast" # Note: blast is in development

# For blast below
path_blast_fasta="" # Ex: /home/shared/edna/reference_database/2024/teleo/blast_db/refdb_all_fish_teleo.fasta

# For ecotag below
path_ecotag_tree="" # Directory not file # Caution: here it needs to be the path to the tree (NCBI tree ; -t ecotag option not -d). Ex: /home/shared/edna/reference_database/2023_06/teleo_custom_embl/customtaxonomy/
path_ecotag_fasta="" # obitools fasta file. Ex: /home/shared/edna/reference_database/2023_06/teleo_custom_embl/db_teleo_custom_and_embl.fasta

# For dada2 below
path_dada_all="" # DADA2 all fasta file. Ex: /home/shared/edna/reference_database/2025/teleo/CRABS_MITOFISH_MIDORI2_GB265_teleo_dada2_all.fasta
path_dada_species="" # DADA2 species fasta file. Ex: /home/shared/edna/reference_database/2025/teleo/CRABS_MITOFISH_MIDORI2_GB265_teleo_dada2_species.fasta

# For decipher below
path_decipher_trained="" # RDS file with trained decipher file. No longer maintained. 

# General
CORES=12