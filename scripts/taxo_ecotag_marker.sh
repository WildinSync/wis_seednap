# Env
# conda activate ref_database

# Save marker name
marker=$3

# Code
ecotag -t $1 -R $2 outputs/02_dada2/$marker/query.fasta > outputs/03_taxo/$marker/query_ecotag.fasta

obiannotate --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
        --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
        --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
        --delete-tag=obiclean_head --delete-tag=obiclean_headcount \
        --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=obiclean_status \
        --delete-tag=seq_length_ori --delete-tag=sminL --delete-tag=sminR \
        --delete-tag=reverse_score --delete-tag=reverse_primer --delete-tag=reverse_match --delete-tag=reverse_tag \
        --delete-tag=forward_tag --delete-tag=forward_score --delete-tag=forward_primer --delete-tag=forward_match \
        --delete-tag=tail_quality outputs/03_taxo/$marker/query_ecotag.fasta > outputs/03_taxo/$marker/query_ecotag_temp.fasta

obitab -o outputs/03_taxo/$marker/query_ecotag_temp.fasta  > outputs/03_taxo/$marker/query_ecotag.tsv

# This outputs the table of sequence - taxonomy correspondance
# Last step is to link each sequence to its abundance table from the dada2 outputs to generate final file

# Env
# conda activate metabarcoding