# Pipeline Steps

## Overview

```
demultiplex (optional) --> trim --> cluster (DADA2 or SWARM) --> taxonomy --> clean (optional) --> export --> report
```

Each step reads outputs from the previous step. The pipeline tracks state in a JSON file, so you can resume from any failed step with `--resume`.

---

## 0. Demultiplex (optional)

**Tool:** Built-in (Cutadapt under the hood)
**Input:** Raw library FASTQ files + a metadata CSV mapping samples to
ligation tags
**Output:** Per-sample FASTQ files under `outputs/01_trim/{marker}/demux/`

The ligation protocol generates per-sample tag files, demultiplexes the
library, detects primers in both orientations, and merges/realigns reads.

**Robustness.** Each sample is processed inside its own
`try`/`except`. A single bad sample does not crash the whole library;
its error is logged and the run continues. If more than
`demultiplex.max_sample_failure_rate` (default 0.5) of samples fail, the
step aborts so an actually-broken library does not silently emit a
mostly-empty output.

**When it runs.** Demultiplexing runs only if `demultiplex` is listed in
`pipeline.steps` (before `trim`). When the raw inputs already arrive as
one FASTQ per sample (most external collaborators), simply omit
`demultiplex` from `steps`. The `standard` protocol is reserved for
future work and currently raises a `NotImplementedError` with a pointer
to the ligation path.

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

### DADA2-by-library

By default (`dada2.per_library: false`) DADA2 learns a single pooled error
model across all input samples -- the unchanged legacy behavior. When
`dada2.per_library: true`, DADA2 instead groups samples by sequencing library
(the FAIRe manifest's `seq_run_id`), learns a separate error model and denoises
each library independently, then merges the per-library sequence tables and
collapses identical ASVs (`mergeSequenceTables` + `collapseNoMismatch`). Use it
for runs that span multiple sequencing runs, where run-specific error profiles
would otherwise be averaged together; it is a no-op for single-library
datasets.

---

## 3. Taxonomic Assignment

**Input:** Representative sequences (`query.fasta`) and abundance table (`otu_table.csv` from SWARM, or `seqtab_clean.csv` from DADA2)
**Output:** Taxonomy CSV in `outputs/03_taxo/{marker}/` and final table at `outputs/{marker}_{method}.csv` (the `{method}` token is `blast` / `ecotag` / `decipher` / `dada2RDP`, e.g. `teleo_dada2RDP.csv` for the DADA2 RDP classifier)

All four methods (BLAST, DADA2 RDP, DECIPHER, ecotag) share a common
post-processor (`seednap.utils.taxonomy.link_taxonomy_with_abundance`),
so the output schema is identical regardless of method: same columns,
same null semantics for missing ranks (cascade null), and the same
`is_contaminant_candidate` column when `taxonomy.contaminants` is set.

The DADA2 RDP and DECIPHER paths take the query FASTA explicitly and
work on either DADA2 ASVs or SWARM OTUs -- they no longer require a
`seqtab_clean.rds` produced by the DADA2 step.

**BLAST LCA algorithm.** The BLAST method resolves multi-hit ambiguity with one
of two header-based, offline LCA resolvers selected by
`taxonomy.blast.lca_algorithm`. The default `cascade` keeps the MEGAN-LR
top-bitscore band (`top_bitscore_pct`, default 10) with a percent-identity floor
(`lca_pident_delta`, default 1) and per-rank identity thresholds
(species 99 / genus 96 / family 90 / order 80 / class 70). The optional
`collapsed_taxonomy` resolver is the eDNAFlow/OceanOmics %identity-window
collapse-to-LCA: hits within `lca_diff` (default 1) percent-identity points of
the best hit above a hard `lca_pid` floor (default 90) collapse to their LCA. It
reads the CRABS lineage from the reference FASTA headers, needs no NCBI
taxids/taxdump, and does not apply cascade's per-rank thresholds.

See [taxonomy-methods.md](taxonomy-methods.md) for detailed method descriptions.

---

## 3b. Decontamination (optional)

**Tool:** Built-in
**Input:** Taxonomy table from step 3 (`{marker}_{method}.csv`)
**Output:** Cleaned table `outputs/{marker}_{method}_cleaned.csv` and a
per-sample `cleaning_report.csv` in the report directory

Runs only when `clean` is listed in `pipeline.steps` (after a feature step, typically
between `taxonomy` and `export`), decontaminating the table against its negative controls. Control
identity (extraction blanks vs PCR blanks, and the extraction batch each sample
belongs to) is derived from the FAIRe manifest, so controls do not need to be
named by convention.

`cleaning.mode` selects the behavior:

- `flag` (default) -- annotate OTUs/ASVs that appear in negative controls
  without changing any counts. Subtraction stays opt-in because it is
  high-consequence.
- `subtract` -- remove control reads from the associated samples: extraction
  blanks clean their own extraction batch, PCR blanks clean the whole dataset.

Export prefers the cleaned table when this step produced one. If the manifest /
control identity is unavailable, the step is skipped with a `[WARN]` and export
falls back to the uncleaned table. The same logic is available standalone on any
abundance table:

```bash
seednap clean {abundance_csv} {field_metadata_csv} {output_csv} [--mode flag|subtract]
```

---

## 4. GBIF Export

**Tool:** Built-in formatter
**Input:** Taxonomy CSV from step 3
**Output:** GBIF-compatible long-format CSV; downstream, a DarwinCore
occurrence CSV via `seednap create-gbif`.

Transforms the wide-format taxonomy table (samples as columns) into GBIF
long format (one row per sample-OTU observation). Adds `rank`
(species/genus/family) and `taxon` (lowest available name) columns.
Zero-count observations are removed. `is_contaminant_candidate` is
carried through and surfaces in the DarwinCore output as
`contamination_flag`.

See [gbif-export.md](gbif-export.md) for the full DarwinCore publishing workflow.

---

## 5. Run Report

**Tool:** Built-in
**Input:** Cutadapt logs, the cluster output (SWARM `otu_table` / DADA2
`track_reads.csv`), and, for the HTML report, the taxonomy table, the SWARM
`otu_table_full.csv`, the run state JSON, and (optionally) dataset metadata.
**Output:** `outputs/04_report/{marker}/read_tracking.{csv,txt}` and
`report.html` (the report directory is configurable via `report.output_dir`).

This step runs when `report` is listed in `pipeline.steps` (it is in the default
steps, so it runs unless you remove it). It always writes the read-tracking table +
step summary; set `report.html_report: false` to skip just the HTML document.

The read-tracking table records per-sample read/sequence counts at each step
(`raw → trimmed → clustered` for SWARM, `raw → trimmed → filtered → denoised →
merged → nonchim` for DADA2) with a `% retained` column, and emits data-loss
warnings against configurable thresholds. Counts that cannot be measured are
recorded as `NA` (never a silent `0`).

The HTML run report is a single self-contained, scientific-paper-style
file: dataset provenance, read-tracking funnel + per-sample retention, a
taxonomy headline (assignment rate per rank, top taxa), feature QC (chimeras,
length), a control/contamination check, the run timeline, and the full console
run log colorized by level (the same palette as the live console).

The read-tracking table is written after the clustering step; the HTML report
after the full run (so taxonomy/provenance are available). Both can be
regenerated from existing outputs with `seednap report MARKER [--html]`.

See [reporting.md](reporting.md) for full details.

---

## Output Directory Structure

```
outputs/
  01_trim/{marker}/              # Trimmed FASTQ files
    demux/                       #   Demultiplexed FASTQ (ligation demux, if enabled)
  02_swarm/{marker}/             # SWARM outputs
    merged/                      #   Merged reads per sample
    dereplicated/                #   Dereplicated per sample
    logs/                        #   Per-step log files
    otu_table.csv                #   Abundance matrix (for taxonomy)
    otu_table_full.csv           #   Full OTU table with metadata
    query.fasta                  #   Representative sequences
  02_dada2/{marker}/             # DADA2 outputs (if used)
  03_taxo/{marker}/              #   BLAST/taxonomy intermediate files
  04_report/{marker}/            # Read-tracking table + HTML report (configurable dir)
    read_tracking.csv            #   Per-sample counts at each step + % retained
    read_tracking.txt            #   Human-readable table
    report.html                  #   Self-contained HTML run report (on by default)
  {marker}_{method}.csv          # Final taxonomy + abundance table (method = blast/ecotag/decipher/dada2RDP)
  {marker}_{method}_cleaned.csv  # Decontaminated table (if "clean" is in pipeline.steps)
  .{marker}_state.json           # Pipeline state (for resume)
```

## State Management and Resume

The pipeline saves state to `.{marker}_state.json` in the output directory after each step. To resume from a failed run:

```bash
seednap run-pipeline config.yaml --resume
```

Completed steps are skipped. Failed steps are retried. Outputs from completed steps are passed forward automatically.
