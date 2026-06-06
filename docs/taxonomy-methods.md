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
      -task {task}
      -outfmt "6 qseqid sseqid pident length mismatch gapopen
               qstart qend sstart send evalue bitscore"
      -perc_identity {perc_identity}
      -qcov_hsp_perc {qcov_hsp_perc}
      -evalue {evalue}
      -max_target_seqs {max_target_seqs}
    ```

    `task` defaults to `megablast` (word size 28), the right call for
    short, high-identity vertebrate amplicons against curated reference
    databases. Switch to `blastn` (word size 11) for divergent references
    where the family/order tier of hits matters.

3. **Phylogeny extraction:** Taxonomy is parsed from reference FASTA headers. Expected format:

    ```
    >ACCESSION\tKingdom;Phylum;Class;Order;Family;Genus;Species
    ```

    The reference parser hard-fails with a descriptive error on malformed
    headers (wrong tab count, wrong semicolon count) so a corrupt DB
    cannot silently produce empty assignments.

4. **Cascade-null per-rank filtering.** The post-processor walks ranks
   from species down to class. When the hit's percent identity is below
   the threshold for a rank, that rank **and every finer rank** are set
   to null. The output therefore never contains orphan ranks like
   `kingdom=Metazoa, phylum=None, class=Mammalia`.

   Per-rank thresholds (YAML defaults follow Pappalardo 2025,
   *Methods in Ecology and Evolution* 16:2380-2394, with rRNA-marker
   tweaks; family raised vs eDNAFlow):

   | Rank | YAML default | CLI shortcut default |
   |---|---|---|
   | `threshold_species` | 99.0 | 98.0 |
   | `threshold_genus` | 96.0 | 96.0 |
   | `threshold_family` | 90.0 | 86.5 |
   | `threshold_order` | 80.0 | (YAML only) |
   | `threshold_class` | 70.0 | (YAML only) |

5. **MEGAN-LR top-bitscore LCA.** The resolver no longer requires exact
   bitscore ties. All hits within `top_bitscore_pct` (default 10%) of
   the best bitscore are pooled, and disagreeing ranks across that pool
   are nulled (Lowest Common Ancestor). Setting `top_bitscore_pct: 0`
   reverts to the old exact-tie behavior.

6. **Output merging.** Taxonomy is **left-joined** onto the OTU/ASV
   abundance table so that OTUs without any BLAST hit surface in the
   final output as `Unassigned` rows rather than being silently dropped.

7. **Contamination flagging.** If `taxonomy.contaminants` is set, every
   row whose `species` matches one of the listed names gets
   `is_contaminant_candidate=True`. Rows are **never** deleted; the flag
   propagates through the GBIF formatter into the DarwinCore output as
   `contamination_flag` for downstream review.

### Reference Database Format

The reference FASTA must have tab-separated headers with semicolon-delimited taxonomy:

```
>KY213962	Metazoa;Chordata;Actinopteri;NA;Centropomidae;Lates;Lates_calcarifer
CACCGCGGTTATACGAGAGGCCCAAGCTGAC...
```

Exactly 7 semicolon-separated ranks are required: kingdom, phylum, class, order, family, genus, species. Use `NA` for unknown ranks.

Databases built with [CRABS](https://github.com/gjeunen/reference_database_creator) (Jeunen et al., 2023) are compatible out of the box. The 2025 CRABS reference DBs write the literal string `NA` where a rank is unknown. SeeDNAP normalizes `NA` (and `""`/`nan`) to a genuine missing rank at the BLAST formatter, in one place, so **neither** LCA resolver treats `NA` as a real taxon -- no over-collapse onto a phantom shared rank, and no literal `NA` leaking into the export. Missing ranks surface as `Unassigned`.

### Configuration

Example with the production cascade defaults made explicit (software
defaults shown in parentheses where they differ):

```yaml
taxonomy:
  method: "blast"
  contaminants:
    - "Homo_sapiens"
    - "Bos_taurus"
  databases:
    blast:
      fasta: "/path/to/reference.fasta"
      perc_identity: 80.0                    # (default: 80.0)
      qcov_hsp_perc: 80.0                    # (default: 80.0)
      evalue: 1.0e-25                        # (default: 1.0e-25)
      max_target_seqs: 5                     # (default: 5)
      task: "megablast"                      # (default: "megablast")
      threshold_species: 99.0                # (default: 99.0)
      threshold_genus: 96.0                  # (default: 96.0)
      threshold_family: 90.0                 # (default: 90.0)
      threshold_order: 80.0                  # (default: 80.0)
      threshold_class: 70.0                  # (default: 70.0)
      top_bitscore_pct: 10.0                 # (default: 10.0)
```

### Optional: collapsed-taxonomy LCA (eDNAFlow/OceanOmics)

The default LCA resolver is `cascade` (steps 4--6 above): the MEGAN-LR
top-bitscore band (`top_bitscore_pct`, default 10) gated by an in-band
identity floor (`lca_pident_delta`, default 1), then per-rank identity
thresholds. An alternative resolver is selectable via
`lca_algorithm: collapsed_taxonomy`, the %identity-window collapse-to-LCA
used by eDNAFlow and OceanOmics.

It works as follows:

1. Discard every hit below `lca_pid` (default 90.0), a hard percent-identity
   floor.
2. Among the surviving hits, take the best percent identity and keep all hits
   within `lca_diff` (default 1.0) identity points of it. This is the
   "identity window".
3. Collapse the lineages of the windowed hits to their Lowest Common
   Ancestor: ranks on which the windowed hits disagree are nulled (and, as
   everywhere in this pipeline, every finer rank cascades to null).

Like `cascade`, it is **header-based**: the lineage comes from the CRABS
reference FASTA headers (Section *Reference Database Format*). It needs **no**
NCBI taxids and no `taxdump`, so it runs fully offline.

How it differs from `cascade`:

- **No per-rank thresholds.** `collapsed_taxonomy` ignores
  `threshold_species`/`threshold_genus`/`threshold_family`/`threshold_order`/
  `threshold_class`; the only identity controls are the `lca_pid` floor and
  the `lca_diff` window. `cascade` keeps its per-rank thresholds (species 99 /
  genus 96 / family 90 / order 80 / class 70).
- **More permissive at low identity.** A 90--96% hit that `cascade` would
  null below the species/genus thresholds can still be reported (down to its
  windowed LCA) once it clears `lca_pid`.
- **More conservative on disagreement.** Because any disagreement *within the
  identity window* collapses to the LCA, a tight cluster of near-equal hits
  spanning two genera resolves only to family, even at high identity, rather
  than picking a rank by per-rank threshold.

Query coverage is enforced separately at the `blastn` step via
`qcov_hsp_perc`, the same as for `cascade`.

`fishbase_tiered` is accepted by the schema but **not implemented** and raises
if selected.

```yaml
taxonomy:
  method: "blast"
  databases:
    blast:
      fasta: "/path/to/reference.fasta"
      lca_algorithm: "collapsed_taxonomy"    # (default: "cascade")
      lca_pid: 90.0                          # (default: 90.0) hard %identity floor
      lca_diff: 1.0                          # (default: 1.0) identity-window width
```

---

## DADA2 RDP Classifier

Uses the naive Bayesian classifier from DADA2 (Wang et al., 2007; Callahan et al., 2016). Requires R and the `dada2` Bioconductor package.

### Algorithm and Bootstrap Threshold

`assignTaxonomy` returns a per-rank bootstrap confidence (0--100). The
post-processor applies a configurable bootstrap threshold (`bootstrap_threshold`,
default **80**, the Wang 2007 recommendation for short rRNA reads):
ranks below the threshold are nulled, and every finer rank cascades to
null in the same way as the BLAST path. The resulting frame matches the
BLAST schema exactly (same column names, same null semantics, contaminant
flag in the same position) so downstream tooling does not branch on
method.

### Configuration

```yaml
taxonomy:
  method: "dada2"
  contaminants:
    - "Homo_sapiens"
  databases:
    dada2:
      all: "/path/to/dada2_all.fasta"
      species: "/path/to/dada2_species.fasta"
      bootstrap_threshold: 80                # (default: 80)
```

The `all` database provides ranks kingdom through genus. The `species` database adds species-level exact matching.

DADA2 RDP works on both DADA2 ASVs and SWARM OTUs; the runner accepts the
query FASTA explicitly and no longer requires a `seqtab_clean.rds` from
the DADA2 step.

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

DECIPHER results are post-processed through the same shared utility as
BLAST and DADA2 RDP, so the output schema and contaminant flag column
are identical across all four methods.

---

## ecotag (OBITools)

Uses the ecotag algorithm from OBITools (Boyer et al., 2016). Requires an NCBI-format taxonomy tree and a reference sequence database.

**Note:** ecotag requires OBITools v1, which has Python 2 dependencies. It
lives in its own conda env. The runner auto-discovers the binary via
`SEEDNAP_OBITOOLS_BIN`, the active `PATH`, or a set of well-known install
paths -- no manual `conda activate obitools` needed when running through
seednap. Setup details: [ecotag-setup.md](ecotag-setup.md).

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
- Huson, D.H. et al. (2018). MEGAN-LR: new algorithms allow accurate binning and easy interactive exploration of metagenomic long reads and contigs. *Biology Direct*, 13, 6.
- Jeunen, G.J. et al. (2023). crabs -- A software program to generate curated reference databases. *Molecular Ecology Resources*, 23, 725-738.
- Murali, A., Bhargava, A. & Wright, E.S. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. *Microbiome*, 6, 140.
- Pappalardo, P. et al. (2025). A field-standard set of identity thresholds for eDNA metabarcoding taxonomic assignment. *Methods in Ecology and Evolution*, 16, 2380-2394.
- Wang, Q. et al. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences. *Applied and Environmental Microbiology*, 73, 5261-5267.
- Whitmore, K. et al. (2023). Sources of contamination in environmental DNA studies. *Nature Ecology and Evolution*, 7, 1-3.
