# Pipeline Steps

## Overview

```
demultiplex (optional) --> trim --> cluster (DADA2 or SWARM) --> taxonomy --> export
```

Each step reads outputs from the previous step. The pipeline tracks state in a JSON file, so you can resume from any failed step with `--resume`.

---

## 1. Primer Trimming

**Tool:** Cutadapt (Martin, 2011)
**Input:** Raw paired-end FASTQ files (R1/R2)
**Output:** Trimmed FASTQ files in `outputs/01_trim/{marker}/`

### Algorithm

SeeDNAP uses a two-pass trimming approach to handle both orientations:

**Pass 1 (5' end):** Trim forward primer from R1 and reverse primer from R2.

```
cutadapt -j {cores} -e {error_rate} -m {min_length} -O {overlap}
  -g {forward_primer} -G {reverse_primer}
  -o {R1_temp} -p {R2_temp}
  {R1_input} {R2_input}
```

**Pass 2 (3' end):** Trim reverse-complement primers from pass 1 output.

```
cutadapt -j {cores} -e {error_rate} -m {min_length} -O {overlap}
  -a {rev_comp_reverse_primer} -A {rev_comp_forward_primer}
  -o {R1_output} -p {R2_output}
  {R1_temp} {R2_temp}
```

Temporary files from pass 1 are deleted after pass 2 completes. If `discard_untrimmed: true`, reads without detected primers are discarded in pass 1.

### File Naming

The trimmer supports both `.R1.fastq` and `_R1.fastq` naming conventions. Output files use the same naming as input.

---

## 2a. SWARM OTU Clustering

**Tools:** VSEARCH (Rognes et al., 2016), SWARM (Mahe et al., 2015)
**Input:** Trimmed FASTQ pairs from `outputs/01_trim/{marker}/`
**Output:** OTU table and representative sequences in `outputs/02_swarm/{marker}/`

### Algorithm

The SWARM workflow proceeds in 7 steps:

**Step 1: Merge paired-end reads** (per sample)

```
vsearch --fastq_mergepairs {R1} --reverse {R2}
  --fastqout {merged}
  --fastq_maxdiffs {maxdiffs}
  --fastq_minovlen {minovlen}
  --fastq_minmergelen {min_seq_length}
  --fastq_maxns 0
```

The `--fastq_maxns 0` flag filters reads containing ambiguous N bases, which SWARM cannot process.

Empty merged files (blanks, negative controls) are detected and skipped with a warning.

**Step 2: Per-sample dereplication**

```
vsearch --fastx_uniques {merged} --fastaout {dereplicated}
  --sizeout --fasta_width 0 --minuniquesize 1 --relabel_sha1
```

SHA1 relabeling ensures identical sequences get the same ID across samples.

Note: For vsearch < 2.28, `--derep_fulllength` is used instead of `--fastx_uniques`.

**Step 3: Global dereplication**

All per-sample FASTA files are concatenated, then globally dereplicated with `--sizein` to sum abundances.

**Step 4: SWARM clustering**

```
swarm {input} -d {d} -t {threads}
  --usearch-abundance
  --internal-structure {struct}
  -s {stats} --seeds {representatives} -o {swarm}
  [--fastidious --boundary {boundary}]
```

**Step 5: Sort representatives by abundance**

```
vsearch --sortbysize {representatives} --output {sorted} --fasta_width 0
```

**Step 6: De novo chimera detection**

```
vsearch --uchime_denovo {sorted} --uchimeout {uchime}
```

**Step 7: Build OTU table**

Parses SWARM cluster membership, per-sample abundances, and chimera status to produce:

- `otu_table_full.csv` -- Complete OTU table with metadata (OTU ID, total reads, chimera status, per-sample counts)
- `otu_table.csv` -- Abundance matrix (sequences x samples), chimeric OTUs removed
- `query.fasta` -- Non-chimeric representative sequences for taxonomy

---

## 2b. DADA2 ASV Processing

**Tool:** DADA2 (Callahan et al., 2016) via R
**Input:** Trimmed FASTQ pairs
**Output:** ASV table and sequences in `outputs/02_dada2/{marker}/`

### Algorithm

1. **Quality filtering:** Filter by expected errors (`maxEE`), truncation quality (`truncQ`), max N bases, PhiX removal
2. **Error learning:** Learn error rates from the data
3. **Denoising:** Infer exact ASVs per sample
4. **Merging:** Merge paired-end reads with minimum overlap requirement
5. **Chimera removal:** De novo chimera detection (consensus or pooled)
6. **ASV table:** Sequence table with per-sample counts

---

## 3. Taxonomic Assignment

**Input:** Representative sequences (`query.fasta`) and abundance table (`otu_table.csv`)
**Output:** Taxonomy CSV in `outputs/03_taxo/{marker}/` and final table at `outputs/{marker}_{method}.csv`

See [taxonomy-methods.md](taxonomy-methods.md) for detailed method descriptions.

---

## 4. GBIF Export

**Tool:** Built-in formatter
**Input:** Taxonomy CSV from step 3
**Output:** GBIF-compatible long-format CSV

Transforms the wide-format taxonomy table (samples as columns) into GBIF long format (one row per sample-OTU observation). Adds `rank` (species/genus/family) and `taxon` (lowest available name) columns. Zero-count observations are removed.

See [gbif-export.md](gbif-export.md) for the full DarwinCore publishing workflow.

---

## Output Directory Structure

```
outputs/
  01_trim/{marker}/              # Trimmed FASTQ files
  02_swarm/{marker}/             # SWARM outputs
    merged/                      #   Merged reads per sample
    dereplicated/                #   Dereplicated per sample
    logs/                        #   Per-step log files
    otu_table.csv                #   Abundance matrix (for taxonomy)
    otu_table_full.csv           #   Full OTU table with metadata
    query.fasta                  #   Representative sequences
  02_dada2/{marker}/             # DADA2 outputs (if used)
  03_taxo/{marker}/              #   BLAST/taxonomy intermediate files
  {marker}_{method}.csv          # Final taxonomy + abundance table
  .{marker}_state.json           # Pipeline state (for resume)
```

## State Management and Resume

The pipeline saves state to `.{marker}_state.json` in the output directory after each step. To resume from a failed run:

```bash
seednap run-pipeline config.yaml --resume
```

Completed steps are skipped. Failed steps are retried. Outputs from completed steps are passed forward automatically.
