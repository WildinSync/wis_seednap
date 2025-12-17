# Save marker name
path_blast_fasta=$1
marker=$2

# Check that blast db files exist - else create them 
if [[ -f "${path_blast_fasta}.njs" ]]; then
    echo "BLAST db files exist."
else
    echo "BLAST db files do not exist. Creating them."
    makeblastdb -dbtype nucl -in ${path_blast_fasta}
fi

# Execute blastn command 
perc_identity=${perc_identity:-80}
qcov_hsp_perc=${qcov_hsp_perc:-80}
evalue=${evalue:-1e-25}
blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"\
    -max_target_seqs 5\
    -perc_identity  "$perc_identity"\
    -qcov_hsp_perc "$qcov_hsp_perc"\
    -evalue "$evalue"\
    -db "${path_blast_fasta}"\
    -query "outputs/02_dada2/${marker}/query.fasta"\
    -out "outputs/03_taxo/${marker}/output_blastn.tsv"

# Execute LCA script from blastn outputs
python scripts/format_blast_outputs.py\
    -b "outputs/03_taxo/${marker}/output_blastn.tsv"\
    -r "${path_blast_fasta}"\
    -c "outputs/02_dada2/${marker}/seqtab_clean.csv"\
    -s "outputs/02_dada2/${marker}/query.fasta"\
    --thre_species 100

# Now rename and move final table
mv outputs/03_taxo/${marker}/output_blastn_final_table.csv outputs/${marker}_blast.csv