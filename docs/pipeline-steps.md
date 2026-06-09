# Pipeline Steps

How each SeeDNAP pipeline step works, which config keys drive it, and what it writes.

The pipeline runs the stages listed in `pipeline.steps`, in order. A stage runs only if it is listed. Valid stages are `demultiplex`, `trim`, `dada2`, `swarm`, `taxonomy`, `clean`, `export`, `report`; `dada2` and `swarm` are mutually exclusive. The default `pipeline.steps` is `[trim, dada2, taxonomy, export, report]`.

```text
demultiplex (optional) --> trim --> dada2 OR swarm --> taxonomy --> clean (optional) --> export --> report
```

Each step reads the previous step's outputs. The pipeline records progress in a state JSON, so a failed run can resume from the failed step with `--resume`.

> [!TIP]
> Run `seednap validate config.yaml` before `run-pipeline`. It runs preflight checks (referenced files exist, taxonomy DB blocks resolve) and catches config errors before any compute starts. `run-pipeline` runs the same preflight automatically. To decode any error code shown by the pipeline, run `seednap explain <code>`.

See [configuration.md](configuration.md) for the full config reference and [cli-reference.md](cli-reference.md) for every command and flag.

## 0. Demultiplex (optional)

Tool: Cutadapt (tag generation + tag matching). Input: one multiplexed library FASTQ pair plus a sample-tag metadata CSV. Output: per-sample FASTQ pairs under `outputs/01_trim/{marker}/demux/`.

The ligation protocol generates per-sample tag files, splits the library by tag, detects primers in both orientations, and realigns reads.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `demultiplex.protocol` | `ligation` \| `standard` \| `none` | `none` | Demultiplexing protocol. Only `ligation` is implemented. |
| `demultiplex.metadata` | path | `None` | CSV mapping samples to ligation tags. Required for the ligation path. |
| `demultiplex.max_sample_failure_rate` | float | `0.5` | Abort the step if more than this fraction of samples fail. |

Each sample is processed in its own `try`/`except`: one bad sample is logged and skipped, not fatal. If more than `demultiplex.max_sample_failure_rate` of samples fail, the step aborts so a broken library does not emit a mostly-empty output.

> [!IMPORTANT]
> The ligation path requires `demultiplex.metadata`. The CSV must contain the columns `eventID`, `tag_demultiplex`, and `library`. A missing `demultiplex.metadata` raises at runtime.

> [!WARNING]
> Listing `demultiplex` in `pipeline.steps` with any protocol other than `ligation` (including the default `none` and the unimplemented `standard`) is rejected at config load, before any step runs. If your reads already arrive as one FASTQ pair per sample (common for external collaborators), omit `demultiplex` from `pipeline.steps` so the pipeline starts at `trim`.

## 1. Primer Trimming

Tool: Cutadapt (Martin, 2011). Input: paired-end FASTQ files (R1/R2). Output: trimmed FASTQ pairs in `outputs/01_trim/{marker}/`.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `trimming.min_length` | int | `20` | Minimum read length after trimming (`-m`). |
| `trimming.max_error_rate` | float | `0.1` | Maximum error rate for primer matching (`-e`). |
| `trimming.overlap` | int | `3` | Minimum overlap for primer detection (`-O`). |
| `trimming.cores` | int | `1` | CPU cores for cutadapt (`-j`). |
| `trimming.discard_untrimmed` | bool | `True` | Discard reads whose 5' primer was not found (pass 1). |

### Algorithm

Two-pass trimming handles both primer orientations.

Pass 1 (5' end): trim forward primer from R1 (`-g`) and reverse primer from R2 (`-G`).

```text
cutadapt -j {cores} -e {max_error_rate} -m {min_length} -O {overlap}
  -g {forward_primer} -G {reverse_primer}
  -o {R1_temp} -p {R2_temp}
  {R1_input} {R2_input}
```

Pass 2 (3' end): trim reverse-complement primers from the pass 1 output (`-a`/`-A`).

```text
cutadapt -j {cores} -e {max_error_rate} -m {min_length} -O {overlap}
  -a {rev_comp_reverse_primer} -A {rev_comp_forward_primer}
  -o {R1_output} -p {R2_output}
  {R1_temp} {R2_temp}
```

Pass 1 temporary files are deleted after pass 2. When `trimming.discard_untrimmed: true`, reads without a detected 5' primer are dropped in pass 1.

### File naming

The trimmer detects inputs in both `.R1.fastq` and `_R1.fastq` conventions (plus `_R1_001.fastq` and `.gz` variants). Trimmed outputs are always written as `{sample}.R1.fastq` / `{sample}.R2.fastq`.

## 2a. SWARM OTU clustering

Tools: VSEARCH (Rognes et al., 2016), SWARM (Mahe et al., 2015). Input: trimmed FASTQ pairs from `outputs/01_trim/{marker}/`. Output: OTU table and representative sequences in `outputs/02_swarm/{marker}/`.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `swarm.merge.fastq_maxdiffs` | int | `10` | Max differences in the overlap region when merging pairs. |
| `swarm.merge.fastq_minovlen` | int | `10` | Min overlap length for merging. |
| `swarm.merge.allow_stagger` | bool | `False` | Allow merging staggered reads (`--fastq_allowmergestagger`). |
| `swarm.clustering.d` | int | `1` | Clustering distance threshold. |
| `swarm.clustering.fastidious` | bool | `True` | Refine singletons (fastidious mode). |
| `swarm.clustering.boundary` | int | `3` | Min mass for large OTUs in fastidious mode. |
| `swarm.clustering.threads` | int | `4` | SWARM threads. |
| `swarm.chimera.method` | `denovo` \| `none` | `denovo` | De novo chimera detection, or skip it. |
| `swarm.min_sequence_length` | int | `20` | Min merged sequence length. |

### Algorithm

The SWARM workflow runs in 7 steps.

Step 1: merge paired-end reads (per sample).

```text
vsearch --fastq_mergepairs {R1} --reverse {R2}
  --fastqout {merged}
  --fastq_maxdiffs {swarm.merge.fastq_maxdiffs}
  --fastq_minovlen {swarm.merge.fastq_minovlen}
  --fastq_minmergelen {swarm.min_sequence_length}
  --fastq_maxns 0
  [--fastq_allowmergestagger]   # when swarm.merge.allow_stagger
```

`--fastq_maxns 0` drops reads with ambiguous N bases, which SWARM cannot process. Empty merged files (blanks, negative controls) are detected and skipped with a warning.

Step 2: per-sample dereplication.

```text
vsearch --fastx_uniques {merged} --fastaout {dereplicated}
  --sizeout --fasta_width 0 --minuniquesize 1 --relabel_sha1
```

SHA1 relabeling gives identical sequences the same ID across samples. For vsearch < 2.28, `--derep_fulllength` is used instead of `--fastx_uniques`.

Step 3: global dereplication. All per-sample FASTAs are concatenated, then globally dereplicated with `--sizein` to sum abundances.

Step 4: SWARM clustering.

```text
swarm {input} -d {swarm.clustering.d} -t {swarm.clustering.threads}
  --usearch-abundance
  --internal-structure {struct}
  -s {stats} --seeds {representatives} -o {swarm}
  [--fastidious --boundary {swarm.clustering.boundary}]   # when fastidious
```

> [!IMPORTANT]
> SWARM's fastidious mode requires `swarm.clustering.d = 1`. With `fastidious: true` and `d > 1`, SWARM exits non-zero and the runner reports a config-mismatch error. Either set `d = 1` or set `fastidious: false`. Run `seednap explain <code>` on the reported error code for the full fix.

Step 5: sort representatives by abundance.

```text
vsearch --sortbysize {representatives} --output {sorted} --fasta_width 0
```

Step 6: de novo chimera detection (only when `swarm.chimera.method: denovo`, the default; `none` skips it).

```text
vsearch --uchime_denovo {sorted} --uchimeout {uchime}
```

Step 7: build the OTU table. Parses cluster membership, per-sample abundances, and chimera status into:

- `otu_table_full.csv`: complete OTU table with metadata (OTU ID, total reads, chimera status, per-sample counts).
- `otu_table.csv`: abundance matrix (sequences x samples), chimeric OTUs removed.
- `query.fasta`: non-chimeric representative sequences for taxonomy.

## 2b. DADA2 ASV processing

Tool: DADA2 (Callahan et al., 2016) via R. Input: trimmed FASTQ pairs. Output: ASV table and sequences in `outputs/02_dada2/{marker}/`.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `dada2.filter.max_ee` | float | `2.0` | Max expected errors (`maxEE`). |
| `dada2.filter.trunc_q` | int | `11` | Truncate reads at the first base with quality <= this (`truncQ`). |
| `dada2.filter.max_n` | int | `0` | Max N bases allowed (`maxN`). |
| `dada2.filter.rm_phix` | bool | `True` | Remove PhiX reads. |
| `dada2.filter.min_len` | int \| null | `None` | Min read length (optional). |
| `dada2.filter.max_len` | int \| null | `None` | Max read length (optional). |
| `dada2.merge.min_overlap` | int | `20` | Min overlap for merging pairs. |
| `dada2.merge.max_mismatch` | int | `0` | Max mismatches in the overlap region. |
| `dada2.chimera.method` | `consensus` \| `pooled` \| `none` | `consensus` | De novo chimera detection mode (or skip). |
| `dada2.pool` | bool | `False` | Pool samples for denoising. |
| `dada2.multithread` | bool | `True` | Use multithreading. |
| `dada2.collect_metrics` | bool | `True` | Write ASV summary stats (DADA2 path only). |
| `dada2.per_library` | bool | `False` | Learn error models per sequencing library, then merge. |

### Algorithm

1. Quality filtering: filter by expected errors (`max_ee`), truncation quality (`trunc_q`), max N bases (`max_n`), optional length bounds, and PhiX removal (`rm_phix`).
2. Error learning: learn error rates from the data.
3. Denoising: infer exact ASVs per sample.
4. Merging: merge pairs with `min_overlap` / `max_mismatch`.
5. Chimera removal: de novo detection per `chimera.method` (`consensus` or `pooled`; `none` skips).
6. ASV table: sequence table with per-sample counts.

When `dada2.collect_metrics: true` (default), ASV summary statistics are written to `outputs/02_dada2/{marker}/metrics/metrics.json` and `metrics.csv`. This is the DADA2 path only; per-step read counts live in the run report, not here.

### DADA2 per-library

- Config key: `dada2.per_library` (default `false`).
- Default (`false`): DADA2 learns one pooled error model across all input samples (the legacy behavior).
- When `true`: DADA2 groups samples by sequencing library (the FAIRe manifest's `seq_run_id`), learns and denoises each library separately, then merges the per-library tables and collapses identical ASVs (`mergeSequenceTables` + `collapseNoMismatch`).
- When to use: runs spanning multiple sequencing runs, where run-specific error profiles would otherwise be averaged together. It is a no-op for single-library datasets.

## 3. Taxonomic assignment

Input: representative sequences (`query.fasta`) and an abundance table (`otu_table.csv` from SWARM, or `seqtab_clean_t.csv` from DADA2). Output: a taxonomy CSV in `outputs/03_taxo/{marker}/` and a final table `outputs/{marker}_{token}.csv`.

The final-table token depends on the method: `blast`, `ecotag`, `decipher`, or `dada2RDP` for the DADA2 RDP classifier (for example `teleo_dada2RDP.csv`).

> [!NOTE]
> The taxonomy table uses the token `dada2RDP` for the DADA2 method, but the cleaned and GBIF tables (sections 3b and 4) use the raw `taxonomy.method` enum value `dada2`. So the DADA2 cleaned table is `{marker}_dada2_cleaned.csv`, not `{marker}_dada2RDP_cleaned.csv`.

All four methods (BLAST, DADA2 RDP, DECIPHER, ecotag) share a post-processor (`seednap.utils.taxonomy.link_taxonomy_with_abundance`), so the output schema is identical regardless of method: same columns, same cascade-null semantics for missing ranks, and the same `is_contaminant_candidate` column when `taxonomy.contaminants` is set. The DADA2 RDP and DECIPHER paths take the query FASTA explicitly and work on either DADA2 ASVs or SWARM OTUs; they do not require a `seqtab_clean.rds`.

### BLAST tuning keys

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `taxonomy.databases.blast.perc_identity` | float | `80.0` | Minimum percent identity passed to blastn. |
| `taxonomy.databases.blast.qcov_hsp_perc` | float | `80.0` | Minimum query coverage per HSP. |
| `taxonomy.databases.blast.evalue` | float | `1e-25` | Maximum e-value. |
| `taxonomy.databases.blast.max_target_seqs` | int | `5` | Maximum target sequences per query. |
| `taxonomy.databases.blast.task` | `megablast` \| `blastn` \| `dc-megablast` \| `blastn-short` | `megablast` | blastn task type. |
| `taxonomy.databases.blast.lca_algorithm` | `cascade` \| `collapsed_taxonomy` \| `fishbase_tiered` | `cascade` | LCA resolver for multi-hit ambiguity. |

### BLAST LCA resolvers

The BLAST method resolves multi-hit ambiguity with one of three header-based, offline LCA resolvers, selected by `taxonomy.databases.blast.lca_algorithm`:

- `cascade` (default): keeps the MEGAN-LR top-bitscore band (`top_bitscore_pct`, default 10) with a percent-identity floor (`lca_pident_delta`, default 1), then applies per-rank identity thresholds (`threshold_species` 99 / `threshold_genus` 96 / `threshold_family` 90 / `threshold_order` 80 / `threshold_class` 70).
- `collapsed_taxonomy`: the eDNAFlow/OceanOmics %identity-window collapse. Hits within `lca_diff` (default 1) percent-identity points of the best hit, above a hard `lca_pid` floor (default 90), collapse to their LCA. It reads the CRABS lineage from the FASTA headers, needs no NCBI taxids/taxdump, and does not apply cascade's per-rank thresholds.
- `fishbase_tiered`: reserved and not yet implemented; selecting it raises `NotImplementedError` at runtime. Use `cascade` or `collapsed_taxonomy`.

See [taxonomy-methods.md](taxonomy-methods.md) for full method descriptions and threshold semantics.

## 3b. Decontamination (optional)

Tool: built-in. Input: the taxonomy table from step 3 (`{marker}_{taxonomy.method}.csv`). Output: a cleaned table `outputs/{marker}_{taxonomy.method}_cleaned.csv` and a per-sample `cleaning_report.csv` in the report directory.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `cleaning.mode` | `flag` \| `subtract` | `flag` | `flag` annotates control OTUs/ASVs without changing counts; `subtract` removes control reads. |

Runs only when `clean` is in `pipeline.steps` (after a feature step, typically between `taxonomy` and `export`). It decontaminates the table against its negative controls. Control identity (extraction blanks vs PCR blanks, and each sample's extraction batch) is derived from the FAIRe manifest, so controls do not need to be named by convention.

> [!NOTE]
> `cleaning.mode` defaults to `flag`, which annotates control features without changing any counts. `subtract` mutates read counts (extraction blanks clean their own extraction batch; PCR blanks clean the whole dataset). Subtraction is high-consequence and stays opt-in. If the manifest or control identity is unavailable, the step is skipped with a `[WARN]` and export falls back to the uncleaned table.

The same logic is available standalone on any abundance table:

```bash
seednap clean {abundance_csv} {field_metadata_csv} {output_csv} \
  [--mode flag|subtract] [--project-metadata PATH] [--id-col COL] [--report PATH]
```

## 4. GBIF export

Tool: built-in formatter. Input: the taxonomy CSV from step 3 (cleaned table preferred when section 3b produced one). Output: a GBIF-compatible long-format CSV `outputs/{marker}_{taxonomy.method}_gbif.csv`; downstream, a DarwinCore occurrence CSV via `seednap create-gbif`.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `export.gbif.add_rank` | bool | `True` | Add a `rank` column (species/genus/family/higher). |
| `export.gbif.add_taxon` | bool | `True` | Add a `taxon` column (lowest available name). |

Transforms the wide-format taxonomy table (samples as columns) into GBIF long format (one row per sample-OTU observation). Zero-count observations are removed. `is_contaminant_candidate` is carried through and surfaces in the DarwinCore output as `contamination_flag`.

See [gbif-export.md](gbif-export.md) for the full DarwinCore publishing workflow.

## 5. Run report

Tool: built-in. Input: Cutadapt logs, the cluster output (SWARM `otu_table` or DADA2 `track_reads.csv`), and, for the HTML report, the taxonomy table, the SWARM `otu_table_full.csv`, the run state JSON, and optional dataset metadata. Output: `read_tracking.{csv,txt}`, `step_summary.csv`, and `report.html` under the report directory (default `outputs/04_report/{marker}/`, configurable via `report.output_dir`).

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `report.html_report` | bool | `True` | Generate the self-contained HTML report. Set `false` to write only the read-tracking table and step summary. |
| `report.output_dir` | path \| null | `None` | Base directory for report artifacts; a per-marker subdir is created inside. Defaults to `<paths.output>/04_report`. |
| `report.warn_below_retention_pct` | float | `30.0` | Warn for samples whose final non-chimeric reads fall below this % of raw reads. |
| `report.warn_step_loss_pct` | float | `70.0` | Warn when a single step drops more than this % of a sample's reads. |
| `report.sample_metadata` | path \| null | `None` | Per-sample field metadata CSV for the report's provenance section. |
| `report.project_metadata` | path \| null | `None` | Project metadata CSV for the report's Dataset section. |

This step runs when `report` is in `pipeline.steps` (it is in the default steps). It always writes the read-tracking table and step summary; `report.html_report: false` skips just the HTML document.

The read-tracking table records per-sample read/sequence counts at each step (`raw -> trimmed -> clustered` for SWARM; `raw -> trimmed -> filtered -> denoised -> merged -> nonchim` for DADA2) with a `% retained` column, and emits data-loss warnings against `report.warn_below_retention_pct` and `report.warn_step_loss_pct`.

> [!NOTE]
> Counts that cannot be measured are recorded as `NA`, never a silent `0`. This is a deliberate correctness guarantee: a silent zero would look like real data loss.

The HTML report is a single self-contained file with:

- dataset provenance,
- the read-tracking funnel and per-sample retention,
- a taxonomy headline (assignment rate per rank, top taxa),
- feature QC (chimeras, length),
- a control/contamination check,
- the run timeline,
- the full console run log, colorized by level.

Both outputs can be regenerated from existing outputs with `seednap report MARKER [--html]`.

See [reporting.md](reporting.md) for full details.

## Output directory structure

```text
outputs/
  01_trim/{marker}/              # Trimmed FASTQ files
    demux/                       #   Demultiplexed FASTQ (ligation demux, if "demultiplex" in pipeline.steps)
  02_swarm/{marker}/             # SWARM outputs
    merged/                      #   Merged reads per sample
    dereplicated/                #   Dereplicated per sample
    logs/                        #   Per-step log files
    otu_table.csv                #   Abundance matrix (for taxonomy)
    otu_table_full.csv           #   Full OTU table with metadata
    query.fasta                  #   Representative sequences
  02_dada2/{marker}/             # DADA2 outputs (if used)
    metrics/                     #   metrics.json / metrics.csv (if dada2.collect_metrics)
  03_taxo/{marker}/              # BLAST/taxonomy intermediate files
  04_report/{marker}/            # Read-tracking + HTML report (dir configurable via report.output_dir)
    read_tracking.csv            #   Per-sample counts at each step + % retained
    read_tracking.txt            #   Human-readable table
    step_summary.csv             #   Run-level reads + feature counts after each step
    report.html                  #   Self-contained HTML run report (when report.html_report)
  {marker}_{token}.csv           # Final taxonomy + abundance (token = blast/ecotag/decipher/dada2RDP)
  {marker}_{taxonomy.method}_cleaned.csv  # Decontaminated table (if "clean" in pipeline.steps); token = dada2/blast/ecotag/decipher
  {marker}_{taxonomy.method}_gbif.csv     # GBIF long table (if "export" in pipeline.steps)
  .{marker}_state.json           # Pipeline state (for resume)
```

> [!NOTE]
> The final taxonomy table uses `dada2RDP` for the DADA2 method, but the `_cleaned.csv` and `_gbif.csv` tables use the raw `taxonomy.method` value `dada2`. See section 3.

## State management and resume

The pipeline saves state to `.{marker}_state.json` in the output directory after each step. To resume a failed run:

```bash
seednap run-pipeline config.yaml --resume
```

Completed steps are skipped, failed steps are retried, and outputs from completed steps are passed forward automatically.

## See also

- [configuration.md](configuration.md) -- full config key reference.
- [cli-reference.md](cli-reference.md) -- every command and flag.
- [taxonomy-methods.md](taxonomy-methods.md) -- taxonomy methods and thresholds.
- [reporting.md](reporting.md) -- run report details.
- [gbif-export.md](gbif-export.md) -- DarwinCore publishing workflow.
