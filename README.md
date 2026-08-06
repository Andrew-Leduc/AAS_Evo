# AAS_Evo

Multi-omics data pipeline for CPTAC3 (Clinical Proteomic Tumor Analysis Consortium). Downloads and processes matched genomics (WXS BAM files from GDC) and proteomics (TMT-labeled RAW files from PDC) data for integrated analysis of amino acid substitutions.

**Current focus:** given a genomic missense mutation, find aberrant amino acid substitutions (AAS, i.e. translation errors) at structurally proximal positions whose combination with the missense is predicted to be strongly interacting — synergistic or aggravating — versus null. Candidate AAS positions come from AlphaFold contact maps; the interaction is scored with **EVE double-mutant predictions** (epistasis = `dmm − (missense_single + swap_single)`). Detected AAS are read out as the within-TMT-set carrier vs non-carrier RAAS delta. Missense stability (AlphaMissense, SPURS ddG) is retained as *annotation*, not a selection filter.

## Repository Structure

```
AAS_Evo/
├── config/
│   └── paths.py                      # Environment-aware path config
├── metadata/                            # Tracked metadata (in-repo)
│   ├── studies.tsv                      # Study registry (tissue → PDC study ID)
│   ├── GDC_meta/                        # GDC manifests & sample metadata
│   │   └── manifests/manifest_wxs_bams.tsv
│   ├── PDC_meta/                        # PDC per-tissue CSVs & consolidated files
│   │   ├── {tissue}/PDC_file_manifest.csv
│   │   ├── {tissue}/PDC_study_biospecimen.csv
│   │   └── {tissue}/PDC_study_experimental.csv
│   ├── Tsour_et_al/                     # External AAS dataset (Tsour et al.)
│   │   ├── pep_to_protein.csv           # SAAP → protein/gene/position mapping
│   │   └── peptide_to_patient.csv       # SAAP → patient/TMT set/RAAS
│   └── mapping_report.tsv              # GDC↔PDC matching report
├── scripts/
│   ├── setup/
│   │   ├── setup_metadata.sh            # Orchestrator: fetch all metadata from APIs
│   │   ├── setup_seq_files.sh           # Download all external reference files
│   │   ├── fetch_gdc_manifest.py        # Query GDC REST API for WXS BAM manifest
│   │   └── fetch_pdc_metadata.py        # Query PDC GraphQL API for per-tissue metadata
│   ├── download/
│   │   ├── gdc/
│   │   │   ├── submit_download.sh       # Submit SLURM download jobs
│   │   │   ├── download_chunk.sh        # Single-chunk SLURM download job
│   │   │   ├── download.py              # GDC download via gdc-client
│   │   │   ├── fetch_metadata.py        # Fetch sample metadata from GDC API
│   │   │   ├── fetch_unmatched_bams.py  # Find WXS BAMs for unmatched PDC patients
│   │   │   └── setup_chunks.sh          # Split manifest into chunks
│   │   ├── pdc/
│   │   │   ├── submit_download.sh       # SLURM job wrapper
│   │   │   ├── download.py              # PDC download (rate-limited)
│   │   │   ├── refresh_urls.py          # Refresh expired PDC signed URLs via API
│   │   │   └── consolidate_metadata.py  # Merge per-tissue CSVs
│   │   └── mapping_report.py            # GDC-PDC sample matching
│   ├── proc_bams/
│   │   ├── run_pipeline.sh              # Wrapper: auto-generates file lists + submits jobs
│   │   ├── submit_variant_call.sh       # SLURM array: variant calling
│   │   ├── submit_vep.sh               # SLURM array: VEP annotation + AlphaMissense
│   │   └── consolidate_missense.sh      # Merge missense mutations
│   ├── mutation_analysis/               # SPURS, AAS target selection, MSA/coevolution
│   │   ├── select_eve_aas_targets.py    # EVE-epistasis AAS target selection (contact-narrowed)
│   │   ├── spurs_predict_per_gene.py    # SPURS ddG stability predictions (+ pLDDT extraction)
│   │   ├── submit_spurs_scores.sh       # SLURM array job for SPURS
│   │   ├── add_spurs_to_missense_table.py  # Merge spurs_ddg into missense table
│   │   ├── backfill_plddt.py            # Add pLDDT to existing ddg_matrix files from cached PDBs
│   │   ├── make_unique_missense_table.py   # Collapse to one row per unique mutation + SPURS + pLDDT
│   │   ├── generate_msas.py
│   │   ├── submit_msa_generation.sh
│   │   ├── build_msa_gene_list.py
│   │   ├── coevolution_analysis.py
│   │   ├── submit_coevolution.sh
│   │   ├── saap_missense_cooccurrence.py   # SAAP × missense co-occurrence (per TMT set)
│   │   ├── saap_cooccurrence_am_analysis.py # AM score distribution for co-occurring mutations
│   │   └── saap_evcoupling_analysis.py     # SAAP × EVCouplings × missense coupling check
│   ├── fasta_gen/                       # Custom proteogenomics FASTAs
│   │   ├── generate_mutant_fastas.py
│   │   ├── combine_plex_fastas.py
│   │   ├── generate_contact_saap_fastas.py  # Per-plex FASTAs of contact-proximal SAAP candidates
│   │   ├── generate_compensatory_fastas.py
│   │   ├── submit_compensatory_fastas.sh
│   │   └── submit_proteogenomics.sh
│   ├── ms_search/                       # FragPipe MS database search
│   │   └── fragpipe/
│   │       ├── run_ms_search.sh         # Orchestrator: setup + submit SLURM array
│   │       ├── generate_manifests.py    # Generate per-plex manifests + patch workflows
│   │       ├── add_decoys.py            # Append rev_ decoy sequences to FASTAs
│   │       ├── submit_fragpipe.sh       # SLURM array job (one plex per task)
│   │       ├── submit_fragpipe_contact.sh   # Same, for the contact-SAAP FASTAs/results
│   │       └── validation/              # Analysis notebooks (contact_saap_analysis.ipynb, ...)
│   └── _retired/                        # Superseded scripts (filter_and_rank.py, maxquant, ...)
└── .claude/
    └── CLAUDE.md                        # Detailed project context
```

## Cluster Data Layout

```
/scratch/leduc.an/AAS_Evo/
├── BAMS/              # GDC BAM files (by UUID subdirectory)
├── RAW/               # PDC RAW files (flattened)
├── VCF/               # Variant calls (per-chunk subdirectories)
│   ├── chunk_00/      # VCFs from first BAM chunk
│   └── chunk_01/ ...
├── VEP/               # VEP annotations (per-chunk subdirectories)
│   ├── chunk_00/      # VEP output from first chunk
│   ├── chunk_01/ ...
│   └── all_missense_mutations.tsv  # Consolidated (top-level)
├── SEQ_FILES/         # Reference files
│   ├── hg38.fa        # Human reference genome (GRCh38)
│   ├── cds.chr.bed    # CDS regions for targeted variant calling
│   ├── uniprot_human_canonical.fasta  # UniProt reviewed proteome
│   ├── uniref50       # MMseqs2 UniRef50 database
│   └── AlphaMissense_hg38.tsv.gz     # AlphaMissense pathogenicity data
├── ANALYSIS/          # Mutation filtering, ranking & SAAP analysis output
│   ├── ranked_mutations.tsv
│   ├── top_5000_mutations.tsv
│   ├── all_missense_with_spurs.tsv        # Full per-sample missense + spurs_ddg
│   ├── unique_missense_mutations.tsv      # One row per unique mutation + spurs_ddg + plddt
│   ├── gene_list_for_msa.txt
│   ├── gene_list_for_spurs_all.txt        # All missense-table genes not yet in SPURS
│   └── saap_cooccurrence/                 # SAAP co-occurrence analysis outputs
│       ├── saap_set_summary.tsv
│       ├── saap_missense_cooccurrence.tsv
│       ├── evc_coverage.tsv
│       ├── evc_cooccurrence.tsv
│       └── evc_summary.txt
├── FASTA/             # Custom proteogenomics FASTAs
│   ├── per_sample/    # Per-sample mutant entries
│   ├── per_plex/      # Reference + plex-specific mutants + compensatory
│   └── compensatory/  # Predicted compensatory mutation entries
├── SPURS/             # SPURS ddG stability predictions
│   ├── pdb_cache/     # AlphaFold PDB files (cached per gene)
│   ├── ddg_matrices/  # {ACC}.{GENE}.ddg_matrix.tsv (pos × 20 AA, with plddt column)
│   └── hf_cache/      # HuggingFace model cache
├── MSA/               # Per-gene multiple sequence alignments (A3M)
├── COEVOL/            # Coevolution analysis output
├── MS_SEARCH/         # FragPipe MS search setup & results
│   ├── manifests/     # Per-plex .fp-manifest files
│   ├── annotations/   # Per-plex TMT channel annotations
│   └── results/       # Per-plex FragPipe output
└── logs/              # SLURM job logs
```

## Pipeline Overview

```
Download BAMs → Variant Call → VEP (+ AlphaMissense) → Consolidate → All Missense Mutations
                                                                            │
              ┌──────────────────────────────────────────────────────────────┤
              │                                                              │
   SPURS ddG (per gene, +pLDDT)                          AlphaFold contact maps
       (annotation only)                                          │
              │                                       Contact-proximal AAS candidates
              │                                                    │
              │                          EVE double-mutant epistasis scoring, two arms:
              │                          synergistic / deleterious / null
              │                          (select_eve_aas_targets.py)
              │                                                    │
              └──────────────────►  Per-plex contact-SAAP FASTAs  ◄┘
                                    (generate_contact_saap_fastas.py)
                                                 │
                                    FragPipe MS Search (per plex)
                                                 │
                       Within-TMT-set carrier vs non-carrier RAAS delta
                                    (contact_saap_analysis.ipynb)
```

**Retained (secondary) branches**, kept but no longer the primary path:
- SAAP co-occurrence vs the Tsour et al. dataset + EVCouplings check (`saap_*` scripts).
- In-house MSA + coevolution compensatory prediction (`generate_msas` → `coevolution_analysis` → compensatory FASTAs) — **superseded as the compensation signal by EVE double-mutant predictions.**
- Composite stability ranking (`filter_and_rank.py`) — moved to `scripts/_retired/`; missense stability is now annotation, not a selection axis.

## Workflows

### 0. Metadata Setup (One-Time)

Fetches all metadata from GDC and PDC APIs programmatically. No manual portal downloads needed.

```bash
bash scripts/setup/setup_metadata.sh
```

This runs six steps automatically:
1. **GDC manifest** — queries GDC REST API for all CPTAC WXS BAM files (projects from `studies.tsv`)
2. **GDC metadata** — enriches manifest with case/sample info
3. **PDC metadata** — queries PDC GraphQL API for file manifests, biospecimen data, and TMT plex layouts (all 12 tissues)
4. **Consolidate** — merges per-tissue PDC CSVs into unified files
5. **Matching** — cross-references GDC↔PDC samples, generates pruned manifests
6. **Recover unmatched** — finds WXS BAMs for PDC patients missing from the initial GDC manifest, merges them in, and re-runs matching

To add a new tissue/study, add a row to `metadata/studies.tsv` and re-run.

**GDC token**: BAM downloads require dbGaP authorization. Download your token from the [GDC portal](https://portal.gdc.cancer.gov) (username → Download Token) and pass the file path to `submit_download.sh`. The token is NOT needed for this metadata setup step.

### 1. Reference File Setup (One-Time)

Downloads and indexes all external reference files (~100 GB total):

```bash
# As SLURM job (recommended for large downloads):
sbatch scripts/setup/setup_seq_files.sh

# Or skip large files (UniRef90, VEP) for initial testing:
bash scripts/setup/setup_seq_files.sh --skip-large
```

Downloads: hg38 genome, GENCODE CDS regions, UniProt canonical proteome, AlphaMissense data, UniRef90 database, and VEP container + cache.

### 2. Data Download (Chunked)

BAMs are downloaded in chunks of ~400 to stay within storage limits:

```bash
# One-time: split manifest into chunks (400 BAMs each, tissue-only)
bash scripts/download/gdc/setup_chunks.sh 400 path/to/manifest_wxs_bams_tissue.tsv

# Per chunk (repeat for each chunk manifest):
bash scripts/download/gdc/submit_download.sh path/to/chunks/chunk_00.tsv
```

**PDC RAW files** (open access):
```bash
sbatch scripts/download/pdc/submit_download.sh
```

### 3. BAM Processing Pipeline

Processes all BAMs currently in `BAMS/`, outputs to per-chunk subdirectories (`VCF/chunk_00/`, `VEP/chunk_00/`):

```bash
# Step 1: Variant calling (finds BAMs automatically, outputs to VCF/chunk_00/)
bash scripts/proc_bams/run_pipeline.sh variant-call chunk_00

# Step 2: VEP annotation with AlphaMissense (reads VCF/chunk_00/, outputs to VEP/chunk_00/)
bash scripts/proc_bams/run_pipeline.sh vep chunk_00

# Step 3: Delete BAMs, download next chunk, repeat with chunk_01, etc.
rm -rf /scratch/leduc.an/AAS_Evo/BAMS/*
```

VCF/VEP outputs persist across chunks in their subdirectories. Both scripts skip already-processed samples.

**Final output** (`all_missense_mutations.tsv`, 18 columns): sample_id, genomic position, consequence, gene symbol, protein change (HGVSp), amino acid swap, gnomADe_AF, AlphaMissense pathogenicity + class, read depths, VAF.

### 4. Consolidate Mutations

```bash
# After ALL chunks are processed: merge per-sample VEP TSVs
bash scripts/proc_bams/consolidate_missense.sh
```

Output: `VEP/all_missense_mutations.tsv` (18 columns). Merge SPURS scores in (§5) to get `ANALYSIS/all_missense_with_spurs.tsv`, the working table for AAS targeting.

> The old composite pathogenicity ranking (`filter_and_rank.py`, top-5000 by AM/gnomAD/recurrence) is retired to `scripts/_retired/` — missense stability is now annotation, not a selection axis. Retrieve it from there if you still want a stability-ranked gene list for the MSA/coevolution branch.

### 5. SPURS Stability Predictions

SPURS predicts per-position ddG (stability change) for all 20 amino acid substitutions using AlphaFold structures. The ddG matrix includes a `plddt` column (AlphaFold per-residue confidence, 0–100) for quality filtering — only trust ddG at positions with pLDDT > 70.

```bash
# One-time setup
bash scripts/mutation_analysis/setup_spurs_env.sh

# Generate gene list of all genes in missense table not yet computed
awk -F'\t' 'NR>1 {print $7}' /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv \
    | sort -u > /tmp/all_missense_genes.txt
ls /scratch/leduc.an/AAS_Evo/SPURS/ddg_matrices/*.ddg_matrix.tsv \
    | xargs -n1 basename | awk -F'.' '{print $2}' | sort -u > /tmp/done_genes.txt
comm -23 /tmp/all_missense_genes.txt /tmp/done_genes.txt \
    > /scratch/leduc.an/AAS_Evo/ANALYSIS/gene_list_for_spurs_all.txt

# Submit in batches of 1000 (Discovery array limit)
NUM=15621  # replace with actual wc -l output
BATCH=1000
for offset in $(seq 0 $BATCH $((NUM - 1))); do
    remaining=$((NUM - offset))
    size=$((remaining < BATCH ? remaining : BATCH))
    GENE_OFFSET=$offset sbatch --array=1-${size}%20 \
        scripts/mutation_analysis/submit_spurs_scores.sh
done

# Backfill pLDDT into existing ddg_matrix files (no model needed, reads cached PDBs)
srun ... python scripts/mutation_analysis/backfill_plddt.py
```

Output: `SPURS/ddg_matrices/{ACC}.{GENE}.ddg_matrix.tsv` — columns: `pos_1based`, `wt_aa`, `plddt`, `to_A` … `to_Y`

After all SPURS jobs complete, build the lightweight unique mutation table:

```bash
python3 scripts/mutation_analysis/make_unique_missense_table.py \
    --missense /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv \
    --spurs-dir /scratch/leduc.an/AAS_Evo/SPURS \
    -o /scratch/leduc.an/AAS_Evo/ANALYSIS/unique_missense_mutations.tsv
```

Output: one row per unique (gene, swap) with `spurs_ddg`, `plddt`, `am_pathogenicity`, `n_samples`, `mean_vaf`.

### 5b. SAAP Co-occurrence Analysis

Cross-references the Tsour et al. AAS dataset with our CPTAC missense calls to ask: for each aberrant amino acid substitution (AAS/SAAP) detected by MS in a TMT set, does any patient in that set carry a genomic missense in the same gene?

**Input data** (`metadata/Tsour_et_al/`):
- `pep_to_protein.csv` — maps each SAAP to its gene, protein, and position
- `peptide_to_patient.csv` — maps each SAAP to the patients/TMT sets it was detected in, with RAAS (relative AAS abundance)

```bash
# Step 1: SAAP × missense co-occurrence (per TMT set)
srun ... python scripts/mutation_analysis/saap_missense_cooccurrence.py \
    --missense /scratch/leduc.an/AAS_Evo/ANALYSIS/all_missense_with_spurs.tsv
# Outputs: saap_set_summary.tsv, saap_missense_cooccurrence.tsv

# Step 2: AM score distribution (co-occurring vs background)
# Co-occurring = missense in SAAP-affected genes, same-set patients
# Background   = missense in non-SAAP genes, same patients & sets
srun ... python scripts/mutation_analysis/saap_cooccurrence_am_analysis.py \
    --missense /scratch/leduc.an/AAS_Evo/ANALYSIS/all_missense_with_spurs.tsv

# Step 3: EVCouplings coupling check
# For co-occurring pairs, checks if the missense position is evolutionarily
# coupled to the SAAP position using precomputed EVChits files (ENSP-level)
srun ... python scripts/mutation_analysis/saap_evcoupling_analysis.py
# Outputs: evc_coverage.tsv, evc_cooccurrence.tsv, evc_summary.txt
```

EVChits files are located at `/projects/marubi/collabs/slavov_rna/evcouplings_files/EVChits/` on Discovery. Each file (`ENSP*_EVChits.csv`) contains evolutionary coupling pairs for SAAP positions in that protein, with `mad_score` and `probability` columns as confidence metrics.

### 5c. Contact-SAAP + EVE-epistasis AAS Targeting (current primary path)

Selects AAS translation-error targets that EVE predicts will interact strongly — or not at all — with an observed missense, then searches for them by MS and reads out the within-TMT-set carrier vs non-carrier RAAS delta.

**Inputs**: `all_missense_with_spurs.tsv`, AlphaFold contact maps (`SPURS/contact_maps/`), and EVE double-mutant matrices at `/projects/marubi/collabs/slavov_rna/aa_mut_preds/eve_alns_double_muts/destab/` (`destab_{ENSP}_{ACC}_b{beta}_double_mutant_matrix.csv`; long format `i, aa_i, j, aa_j, i_score, j_score, dmm_score`, more negative = more deleterious).

```bash
# Select AAS target arms (contacts narrow candidates, EVE scores epistasis)
#   epistasis = dmm_score - (missense_single + swap_single)
#   arms: 2500 synergistic (+ep) / 2500 deleterious (-ep) / 5000 null (|ep|~0)
#   eligibility: missense with >=2 carrier & >=2 non-carrier channels in >=1 TMT set
python3 scripts/mutation_analysis/select_eve_aas_targets.py
# Outputs (ANALYSIS/eve_aas_targets/):
#   eve_aas_targets.tsv          arm-labeled (missense, AAS) pairs  [Stage 1, precomputed EVE]
#   eve_aas_candidates.tsv       all scored candidates
#   eve_proteins_to_compute.tsv  proteins with eligible missense + contacts but no EVE matrix [Stage 2]
#   selection_summary.txt
```

EVE double-mutant matrices are **precomputed for only ~84 proteins**, so Stage 1 arms come mostly from a few dozen genes; `eve_proteins_to_compute.tsv` is the run-EVE list to expand coverage (Stage 2 — re-running the selector folds new matrices in automatically).

The per-plex contact-SAAP FASTAs and their FragPipe search are generated by `generate_contact_saap_fastas.py` and `submit_fragpipe_contact.sh` (results in `MS_SEARCH/results_contact/`). The readout — carrier vs non-carrier RAAS delta computed **within each TMT set**, versus EVE epistasis — lives in `scripts/ms_search/fragpipe/validation/contact_saap_analysis.ipynb`.

> **Gotcha**: the carrier vs non-carrier RAAS delta must be computed within a plex first, then averaged across plexes. Pooling channels across TMT sets mixes batches and *inverts* the sign of the EVE-epistasis correlation.

### 6. MSA Generation (secondary — in-house coevolution branch)

```bash
# Submit MSA generation (auto-finds gene list from ANALYSIS/)
NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/ANALYSIS/gene_list_for_msa.txt)
sbatch --array=1-${NUM_GENES}%10 scripts/mutation_analysis/submit_msa_generation.sh
```

MSA files named by UniProt accession (`P04637.a3m`). Pre-existing MSAs are auto-detected and skipped.

### 7. Coevolution Analysis

Predicts compensatory translation errors using EVcouplings/Direct Coupling Analysis (DCA). Given a destabilizing missense mutation at position i, uses the Potts model coupling tensor J(i,a;j,b) to identify covarying positions and predict which amino acid substitution could compensate.

**Method**: Mean-field DCA computes evolutionary couplings between all position pairs. The coupling parameters directly encode which amino acid combinations are evolutionarily preferred, enabling more accurate compensatory predictions than simple mutual information.

```bash
# After MSA generation completes (auto-finds gene list from ANALYSIS/):
sbatch scripts/mutation_analysis/submit_coevolution.sh
```

Output: `COEVOL/compensatory_predictions.tsv`

### 8. Generate Compensatory FASTAs

```bash
# After coevolution analysis completes:
sbatch scripts/fasta_gen/submit_compensatory_fastas.sh
```

Output: `FASTA/compensatory/all_compensatory.fasta` — tryptic peptides containing both the original destabilizing mutation and the predicted compensatory substitution.

### 9. Proteogenomics FASTA Generation

```bash
# Generate per-sample mutant FASTAs + per-plex search databases
# Automatically includes compensatory entries if FASTA/compensatory/ exists
sbatch scripts/fasta_gen/submit_proteogenomics.sh
```

Creates custom MS search databases per TMT plex: reference proteome + plex-specific mutant tryptic peptides + plex-specific compensatory peptides.

**Tryptic peptide approach**: Instead of adding full mutant proteins (~500 AA each), the pipeline extracts only tryptic peptides containing mutations (~15 AA avg). This minimizes database size and improves FDR statistics.

**FASTA headers** — two formats exist at different stages:

Per-sample FASTAs (`FASTA/per_sample/`) use a simple internal format:
```
>mut|P04637|TP53|R273H|genetic
SVTCTYSPALNKMFCQLAK
```

Per-plex FASTAs (`FASTA/per_plex/`) are rebuilt by `combine_plex_fastas.py` into a Philosopher-compatible mock-UniProt format required for FragPipe:
```
>sp|P04637-R273H-A3F2|TP53-mut TP53 mutant R273H OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=1
>sp|P04637-comp-R273H-G245S-A3F2|TP53-comp TP53 compensatory R273H_G245S OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=1
```

The accession field structure is `{UniProtID}-{swap}-{4-char-seq-hash}`. The hash disambiguates multiple tryptic peptides for the same mutation (different missed cleavage contexts). See the FASTA structure section below for why this format is required.

### 10. MS Database Search (FragPipe)

FragPipe is the primary tool for the MS database search. Each TMT plex is searched independently against its own custom FASTA (reference + plex-specific mutant peptides + compensatory peptides).

#### FragPipe Installation

FragPipe requires several components to be installed together. The pipeline was developed and validated with **FragPipe 24.0**.

Required components (all managed via the FragPipe GUI on first install):
- **FragPipe 24.0** — main orchestrator (`fragpipe-24.0/bin/fragpipe`)
- **MSFragger** — database search engine (bundled with FragPipe)
- **Philosopher** — statistical validation and reporting (bundled with FragPipe)
- **IonQuant / TMT-Integrator** — TMT reporter ion quantification (bundled with FragPipe)
- **Java 17** — required runtime (`jdk-17.0.18+8` used here; set `JAVA_HOME`)

> **Note**: The exact versions of MSFragger, Philosopher, and IonQuant bundled with FragPipe 24.0 should be noted here — check the FragPipe GUI's Tools tab for what was actually used.

Installation paths are configured in `submit_fragpipe.sh`:
```bash
FRAGPIPE_BIN=/home/leduc.an/bin/fragpipe-24.0/bin/fragpipe
FRAGPIPE_TOOLS=/home/leduc.an/bin/fragpipe-24.0/tools
JAVA_HOME=/home/leduc.an/bin/jdk-17.0.18+8
```

#### FASTA Database Structure

FragPipe requires decoy sequences in the search database (target-decoy FDR estimation). The pipeline maintains two separate copies of the per-plex FASTAs:

- `FASTA/per_plex/` — clean FASTAs, no decoys (used by MaxQuant if needed)
- `FASTA/per_plex_fragpipe/` — FASTAs with reversed decoy entries appended (`rev_` prefix)

`run_ms_search.sh` runs `add_decoys.py` automatically to generate the FragPipe copies. Decoys are full sequence reversals (not shuffled), prefixed with `rev_`. FragPipe detects decoys by this prefix.

The per-plex FASTA contains four entry types:
```
# Standard reference proteins (passed through from UniProt)
>sp|P04637|P53_HUMAN Cellular tumor antigen p53 ... GN=TP53 PE=1 SV=4
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP...

# Observed mutant tryptic peptides (Philosopher-compatible mock-UniProt format)
>sp|P04637-R273H-A3F2|TP53-mut TP53 mutant R273H OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=1
SVTCTYSPALNKMFCQLAK

# Compensatory mutation peptides
>sp|P04637-comp-R273H-G245S-A3F2|TP53-comp TP53 compensatory R273H_G245S OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=1
VVRCPHHERCSDSDGLAPPQHLIRVEGNLHAEYLDDQQFIFHSVVVPYEPPEVGSDCTTIHYNYMCNS

# Reversed decoys (appended by add_decoys.py)
>rev_sp|P04637|P53_HUMAN ...
```

**Why this specific header format:** Philosopher classifies FASTA entries by accession format during `philosopher database --annotate`. Entries with non-UniProt accessions (e.g. `MUT_P04637_R273H`) are classified as "generic" proteins and skipped — they never appear in `db.bin`, causing the filter step to crash with "protein not in database" errors. Using the real UniProt accession as the base (`P04637-R273H-HASH`) causes Philosopher to classify the entry as a UniProt variant, store it in `db.bin`, and handle it correctly throughout the pipeline. The `GN=` field is separately required by TMT-Integrator for gene-level intensity reports. Underscores in the accession field also break Philosopher's parser, which is why compensatory swap combos use dashes (`R273H-G245S`) not underscores in the accession field.

Mutant and compensatory entries are **single tryptic peptides** (~15 AA avg), not full proteins. Any mutant peptide that is a substring of its parent reference protein is filtered out at this stage — MSFragger matches by sequence, so a shared substring would assign the PSM to the reference entry and cause Philosopher to crash looking up the mutant accession.

#### Workflow Configuration

The FragPipe workflow is a `.workflow` file exported from the FragPipe GUI. Key settings used:

- **Search type**: TMT-11 closed search
- **TMT channels**: 11-plex (126–131C); forced via `tmtintegrator.channel_num=TMT-11` patch in `generate_manifests.py`
- **Decoy prefix**: `rev_` (must match what `add_decoys.py` uses)
- **Protein FDR**: default template value (standard `--prot 0.05`)

> **Note**: The specific search parameters in the workflow (precursor/fragment mass tolerances, variable modifications, enzyme, missed cleavages, etc.) are set in the `.workflow` file on the cluster. These should be documented here once confirmed — export the workflow from the FragPipe GUI and check the MSFragger section.

To create the workflow:
1. Open FragPipe GUI on a machine with the RAW files accessible (or use a test file)
2. Load files → select TMT-11 experiment type
3. Configure MSFragger search parameters for closed search
4. Enable Philosopher (PeptideProphet, ProteinProphet) and TMT-Integrator
5. **Export workflow** → save as `fragpipe.workflow`
6. Copy to the cluster and pass to `run_ms_search.sh`

#### Running the Search

```bash
# Step 1: Set up per-plex manifests + copy FASTAs with decoys + submit jobs
bash scripts/ms_search/fragpipe/run_ms_search.sh /path/to/template.workflow

# Or: generate manifests first, add workflow manually, then submit
bash scripts/ms_search/fragpipe/run_ms_search.sh
# -> place workflow at MS_SEARCH/fragpipe.workflow, then:
NUM_PLEXES=$(wc -l < /scratch/leduc.an/AAS_Evo/MS_SEARCH/plex_list.txt)
sbatch --array=1-${NUM_PLEXES}%5 scripts/ms_search/fragpipe/submit_fragpipe.sh
```

`run_ms_search.sh` does three things automatically:
1. Copies per-plex FASTAs to `FASTA/per_plex_fragpipe/` with `rev_` decoys appended
2. Generates per-plex `.fp-manifest` files and TMT channel annotation files
3. Patches per-plex workflow files with the correct FASTA path and TMT-11 channel count

Each SLURM job: 16 CPUs, 64 GB RAM, 24h time limit, 5 plexes concurrent.

#### Output

Results are written to `MS_SEARCH/results/{plex_id}/`. Key output files:
- `{plex_id}_1/psm.tsv` — PSM-level results with TMT reporter intensities
- `combined_protein.tsv` — protein-level summary (Philosopher)
- `experiment_annotation.tsv` — TMT channel → patient ID mapping (copied from `annotations/`)

Search completion is detected by the presence of `combined_protein.tsv`. Re-submitting the array job safely skips completed plexes.

## Reference Files

All reference files are downloaded and indexed by `scripts/setup/setup_seq_files.sh`. They live in `/scratch/leduc.an/AAS_Evo/SEQ_FILES/`:

| File | Source | Size |
|------|--------|------|
| `hg38.fa` + `.fai` | UCSC Genome Browser (chr-prefix, matching GDC BAMs) | ~3 GB |
| `cds.chr.bed` | GENCODE v46 CDS regions (merged, standard chromosomes) | ~2 MB |
| `uniprot_human_canonical.fasta` | UniProt reference proteome UP000005640 (reviewed, canonical) | ~25 MB |
| `AlphaMissense_hg38.tsv.gz` + `.tbi` | DeepMind AlphaMissense pathogenicity predictions | ~6 GB |
| `uniref50` (MMseqs2 db) | UniProt Reference Clusters at 50% identity | ~60 GB |
| VEP container + cache | Ensembl VEP Apptainer image | ~15 GB |

## Requirements

**Genomics pipeline:**
- Python 3.8+, pandas, numpy
- `gdc-client` (GDC Data Transfer Tool, for controlled-access BAM downloads)
- `samtools`, `bcftools` (variant calling)
- `bedtools` (for generating CDS BED from GENCODE annotation)
- `htslib` / `tabix` (AlphaMissense indexing)
- Ensembl VEP via Apptainer container (annotation + AlphaMissense plugin)

**Coevolution pipeline:**
- `mmseqs2` (MSA generation against UniRef50)

**MS search pipeline:**
- FragPipe 24.0 with bundled MSFragger, Philosopher, and IonQuant/TMT-Integrator
- Java 17 (required by FragPipe/MSFragger)
- Workflow file configured for TMT-11 closed search (exported from FragPipe GUI)

## Data Sources

| Source | Data Type | Access | Files |
|--------|-----------|--------|-------|
| [GDC](https://portal.gdc.cancer.gov) | WXS BAM files | Controlled (dbGaP) | ~2,098 matched |
| [PDC](https://pdc.cancer.gov) | TMT RAW files | Open (signed URLs) | ~5,020 matched |

## License

MIT License
