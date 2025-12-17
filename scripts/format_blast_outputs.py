# Author: Théophile Sanchez
# Date: 2024-08-15
# Description: format_blast_outputs.py is used to process outputs from blast by retrieving the phylogeny from the
# reference database. The phylogeny of hits with identical bitscores are merged to the most recent common clade.
# Thresholds can be set up to filter hits with low identity percentages. This scripts need to be adapted depending
# on the format of the reference database.

import pandas as pd
import numpy as np
import argparse


def format_csv(tsv_file, fasta_file, output_file):
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
        'bitscore', 'qseq', 'sseq'
    ]
    tsv_df = pd.read_csv(tsv_file, sep='\t', header=None, names=columns)
    with open(fasta_file, 'r') as f:
        fasta_lines = f.readlines()
    fasta_dict = {}
    for line in fasta_lines:
        if line.startswith('>'):
            key = line.split()[0].strip('>')
            fasta_dict[key] = line

    taxonomic_keys = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for taxo_key in taxonomic_keys:
        tsv_df[taxo_key] = None

    for i, row in tsv_df.iterrows():
        phylo = dict(zip(taxonomic_keys, fasta_dict[row.sseqid].replace('\n', '').split('\t')[1].split(';')))
        for key in taxonomic_keys:
            tsv_df.at[i, key] = phylo[key]
    tsv_df['blast_rank'] = tsv_df.groupby('qseqid').cumcount() + 1
    tsv_df.to_csv(output_file, sep='\t', header=True, index=False)


def filter_phylo(data, thresholds):
    for phylo_level, threshold in thresholds.items():
        data.loc[pd.to_numeric(data['pident']) < float(threshold), phylo_level] = None
    return data


def check_ambiguous_hits(taxonomic_assignment_table):
    taxonomic_assignment_table = taxonomic_assignment_table.reset_index(drop=False)
    if len(taxonomic_assignment_table) > 1:
        phylo = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        best_bit_score = taxonomic_assignment_table['bitscore'].iloc[0]
        ambiguous_hits = taxonomic_assignment_table[taxonomic_assignment_table['bitscore'] == best_bit_score].copy()
        if len(ambiguous_hits) > 1:
            same_phylo = np.all([ambiguous_hits[col][ambiguous_hits[col].notna()].nunique() < 2 for col in phylo])

            if not same_phylo:
                taxonomic_assignment_table['keep_for_analysis'] = False
                combined_row = pd.DataFrame([[
                    ambiguous_hits[col].iloc[0] if ambiguous_hits[col].nunique() == 1 else None
                    for col in ambiguous_hits.columns
                ]],
                                            columns=ambiguous_hits.columns)

                combined_row['keep_for_analysis'] = True
                taxonomic_assignment_table = pd.concat([taxonomic_assignment_table, combined_row], ignore_index=True)

    return taxonomic_assignment_table


def fasta_to_dataframe(fasta_file):
    headers = []
    sequences = []
    with open(fasta_file, 'r') as file:
        sequence = ''
        header = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence:
                    headers.append(header)
                    sequences.append(sequence)
                    sequence = ''
                header = line[1:]  # Remove the '>'
            else:
                sequence += line

        if header and sequence:
            headers.append(header)
            sequences.append(sequence)

    df = pd.DataFrame({'ASV_ID': headers, 'Sequence': sequences})
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--blast_outputs', type=str)
    parser.add_argument('-r', '--ref_db', type=str)
    parser.add_argument('-c', '--asv_count', type=str)
    parser.add_argument('-s', '--asv_seq', type=str)
    parser.add_argument('--thre_species', type=float, default=98)
    parser.add_argument('--thre_genus', type=float, default=96)
    parser.add_argument('--thre_family', type=float, default=86.5)
    args = parser.parse_args()

    asv_formatted_path = f'{args.blast_outputs.split(".")[0]}_formatted.tsv'
    format_csv(args.blast_outputs, args.ref_db, asv_formatted_path)
    asv_count = pd.read_csv(args.asv_count, sep=',', index_col=0).T
    asv_formated = pd.read_csv(asv_formatted_path, sep='\t')
    all_asv_seq = fasta_to_dataframe(args.asv_seq)
    asv_count = pd.merge(asv_count, all_asv_seq, how='inner', left_index=True,
                         right_on='Sequence').drop(columns='Sequence')
    asv_formated['keep_for_analysis'] = asv_formated['blast_rank'] == 1
    phylogenetic_thresholds = {'species': args.thre_species, 'genus': args.thre_genus, 'family': args.thre_family}
    asv_formated = asv_formated.replace('None', None)
    asv_formated = filter_phylo(asv_formated, phylogenetic_thresholds)
    asv_formated = asv_formated.set_index('qseqid').groupby('qseqid').apply(check_ambiguous_hits)
    asv_formated = asv_formated[asv_formated['keep_for_analysis'] == True]
    phylo = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    asv_formated = asv_formated[['qseqid', 'pident'] + phylo].rename(columns={'qseqid': 'ASV_ID'})
    final_table = pd.merge(asv_formated, asv_count, how='inner', on='ASV_ID')
    final_table['asv_num'] = final_table['ASV_ID'].str.extract(r'(\d+)').astype(int)
    final_table = final_table.sort_values('asv_num')
    final_table = final_table.drop(columns='asv_num')
    final_table = pd.merge(final_table, all_asv_seq, how='inner', on='ASV_ID')
    final_table.to_csv(f'{args.blast_outputs.split(".")[0]}_final_table.csv', header=True, index=False)


if __name__ == "__main__":
    main()
