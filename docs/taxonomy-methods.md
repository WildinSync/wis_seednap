# Taxonomic Assignment Methods

SeeDNAP supports four taxonomic assignment methods. Each method is selected via the `taxonomy.method` field in the YAML config.

## Method Comparison

| Method | Algorithm | Speed | Best For |
|---|---|---|---|
| **BLAST + LCA** | Local alignment + Lowest Common Ancestor | Moderate | Custom databases, configurable thresholds |
| **DADA2 RDP** | Naive Bayesian classifier | Fast | Standard workflows, DADA2 format databases |
| **DECIPHER** | IdTaxa machine learning classifier | Fast | Pre-trained models, confidence scores |
| **ecotag** | OBITools global alignment | Slow | Legacy OBITools workflows |

---

## BLAST + LCA (Recommended)

### Algorithm

1. **Database creation:** `makeblastdb` creates a nucleotide BLAST database from the reference FASTA (if not already present).

2. **BLASTn search:** Each query sequence is searched against the database.

    ```
    blastn -query {query} -db {db}
      -outfmt "6 qseqid sseqid pident length mismatch gapopen
               qstart qend sstart send evalue bitscore qseq sseq"
      -perc_identity {perc_identity}
      -qcov_hsp_perc {qcov_hsp_perc}
      -evalue {evalue}
      -max_target_seqs {max_target_seqs}
    ```

3. **Phylogeny extraction:** Taxonomy is parsed from reference FASTA headers. Expected format:

    ```
    >ACCESSION\tKingdom;Phylum;Class;Order;Family;Genus;Species
    ```

4. **Phylogenetic filtering:** Taxonomic ranks are set to null if percent identity is below threshold:
    - Species: requires >= `threshold_species`
    - Genus: requires >= `threshold_genus`
    - Family: requires >= `threshold_family`

5. **LCA resolution:** When multiple hits share the same best bitscore:
    - If all hits agree on taxonomy at all ranks: keep first hit
    - If hits disagree: create a consensus row where disagreeing ranks are set to null (Lowest Common Ancestor)
    - Only the resolved best hit is retained per query

6. **Output merging:** Taxonomy is joined with the OTU/ASV abundance table and representative sequences.

### Reference Database Format

The reference FASTA must have tab-separated headers with semicolon-delimited taxonomy:

```
>KY213962	Metazoa;Chordata;Actinopteri;NA;Centropomidae;Lates;Lates_calcarifer
CACCGCGGTTATACGAGAGGCCCAAGCTGAC...
```

Exactly 7 semicolon-separated ranks are required: kingdom, phylum, class, order, family, genus, species. Use `NA` for unknown ranks.

Databases built with [CRABS](https://github.com/gjeunen/reference_database_creator) (Jeunen et al., 2023) are compatible out of the box.

### Configuration

Example with recommended values for eDNA metabarcoding (software defaults in parentheses where different):

```yaml
taxonomy:
  method: "blast"
  databases:
    blast:
      fasta: "/path/to/reference.fasta"
      perc_identity: 80.0                    # (default: 80.0)
      qcov_hsp_perc: 80.0                   # (default: 80.0)
      evalue: 1.0e-10                        # (default: 1.0e-25)
      max_target_seqs: 10                    # (default: 5)
      threshold_species: 100.0               # (default: 98.0)
      threshold_genus: 96.0                  # (default: 96.0)
      threshold_family: 86.5                 # (default: 86.5)
```

---

## DADA2 RDP Classifier

Uses the naive Bayesian classifier from DADA2 (Wang et al., 2007; Callahan et al., 2016). Requires R and the `dada2` Bioconductor package.

### Configuration

```yaml
taxonomy:
  method: "dada2"
  databases:
    dada2:
      all: "/path/to/dada2_all.fasta"
      species: "/path/to/dada2_species.fasta"
```

The `all` database provides ranks kingdom through genus. The `species` database adds species-level exact matching.

---

## DECIPHER IdTaxa

Uses the DECIPHER IdTaxa classifier (Murali et al., 2018). Requires a pre-trained `.rds` classifier file and the R `DECIPHER` package.

### Configuration

```yaml
taxonomy:
  method: "decipher"
  databases:
    decipher:
      trained: "/path/to/trained_classifier.rds"
      threshold: 60
      processors: 8
```

The `threshold` parameter (0-100) controls confidence required for assignment. Lower values assign more sequences but with less certainty.

---

## ecotag (OBITools)

Uses the ecotag algorithm from OBITools (Boyer et al., 2016). Requires an NCBI-format taxonomy tree and a reference sequence database.

**Note:** ecotag requires OBITools v1, which has Python 2 dependencies. Use a separate conda environment.

### Configuration

```yaml
taxonomy:
  method: "ecotag"
  databases:
    ecotag:
      tree: "/path/to/ncbi/taxonomy/"
      fasta: "/path/to/reference.fasta"
```

---

## References

- Boyer, F. et al. (2016). obitools: a unix-inspired software package for DNA metabarcoding. *Molecular Ecology Resources*, 16, 176-182.
- Callahan, B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581-583.
- Camacho, C. et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421.
- Jeunen, G.J. et al. (2023). crabs -- A software program to generate curated reference databases. *Molecular Ecology Resources*, 23, 725-738.
- Murali, A., Bhargava, A. & Wright, E.S. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. *Microbiome*, 6, 140.
- Wang, Q. et al. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences. *Applied and Environmental Microbiology*, 73, 5261-5267.
