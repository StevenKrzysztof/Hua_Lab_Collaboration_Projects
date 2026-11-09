# Tumor RNA-seq: NT vs sg7 tumors, PBS vs anti-PD1

Analysis code for the bulk tumor RNA-seq in [manuscript title]. 12 tumors: genotype NT or sg7 × treatment PBS or
anti-PD1, n = 3 per group. Main question: the PD1 response within each genotype, and whether sg7 changes it (interaction).

## Repository layout
```
├── README.md
├── metadata/
│   ├── sample_metadata.csv      # sample, genotype, treatment, raw FASTQ folder, sequencing batch, notes
│   └── gene_sets.csv            # custom immune modules + genes of interest (set, gene)
├── config/
│   ├── reference_genome.txt     # genome, annotation, indexes, ID / symbol handling
│   ├── parameters.txt           # every command-line parameter and threshold
│   └── software_versions.txt    # written by scripts/collect_versions.sh
├── scripts/
│   ├── upstream/                # SLURM jobs (cluster)
│   │   ├── 00_build_refs_vM25.sh               # STAR + Salmon indexes, tx2gene (GENCODE vM25)
│   │   ├── 01_STAR_Salmon_vM25.sh              # trim -> validate -> STAR -> Salmon, 10 tumors
│   │   └── 01b_STAR_Salmon_vM25_ND10R2_ND11L2.sh   # same, the 2 tumors from the second delivery
│   ├── R/
│   │   ├── 00_config.R          # paths, samples, thresholds, colours, helpers (sourced by all)
│   │   ├── 03_tximport.R        # Salmon -> gene counts
│   │   ├── 04_DESeq2.R          # model, QC, contrasts, NT-vs-sg7 comparison
│   │   ├── 05_enrichment_ORA_GSEA.R
│   │   ├── 06_figures.R         # volcanoes, heatmaps, module heatmaps + scores
│   │   └── 07_supplementary_tables.R
│   ├── TE/                      # transposable-element analysis (same BAMs, same design)
│   │   ├── 01_TElocal_vM25.sh  02_TEcount_vM25.sh  03_TE_DESeq2_PD1vsPBS.Rmd
│   ├── run_R_pipeline.sh        # runs R steps 03 -> 07 in order
│   └── collect_versions.sh      # records tool + package versions into config/
└── results/                     # run logs only (outputs are written to OUTDIR, see 00_config.R)
```

## Order of execution
| Step | Command | Input | Main output |
|---|---|---|---|
| 0 | `sbatch scripts/upstream/00_build_refs_vM25.sh` | GRCm38 FASTA, vM25 GTF | STAR / Salmon indexes, `tx2gene_vM25.tsv` |
| 1 | `sbatch scripts/upstream/01_STAR_Salmon_vM25.sh` and `01b_...` | raw FASTQ (`metadata/sample_metadata.csv`) | `STAR_vM25/`, `Salmon_Results_vM25/` |
| 2 | `bash scripts/run_R_pipeline.sh` | `Salmon_Results_vM25/*/quant.sf` | everything below |
| TE | `sbatch scripts/TE/01_TElocal_vM25.sh`, `02_TEcount_vM25.sh`, then knit `03_TE_DESeq2_PD1vsPBS.Rmd` | STAR BAMs | TE count tables, TE DE results |
| – | `bash scripts/collect_versions.sh` | – | `config/software_versions.txt` |

R steps can also be run one at a time (`Rscript scripts/R/04_DESeq2.R` from the repo root); each step saves an `.RData`
file in `OUTDIR/rds/` that the next step loads. Paths are set at the top of `scripts/R/00_config.R`
(or with the environment variables `TUMOR_RNASEQ_BASE`, `TUMOR_RNASEQ_OUT`, `TX2GENE_VM25`).

## Outputs (in `OUTDIR`) and the figures they support
| Folder / file | Content |
|---|---|
| `QC/` | PCA (PC1–PC3), sample-distance heatmap |
| `DE/DE_<contrast>_all.csv`, `_sig.csv`, `DE_summary.csv` | DESeq2 results: NT PD1 vs PBS, sg7 PD1 vs PBS, interaction |
| `DE/Compare_NT_vs_sg7_PD1_response.csv`, `plots/Venn_*`, `plots/Quadrant_*` | shared / NT-only / sg7-only PD1 genes |
| `plots/Volcano_*` | volcano plots (NT and sg7 on identical axes) |
| `plots/Heatmap_*` | VST z-score heatmaps (top 50, shared / NT-only / sg7-only, interaction, immune modules) |
| `plots/immune_module_scores.*` | per-sample module scores (custom signatures) |
| `enrichment/ORA/`, `enrichment/GSEA/` | ORA and GSEA tables + plots; `GSEA_NES_compare_*` = NT vs sg7 vs interaction |
| `supplementary/` | raw counts, TPM, normalized counts, VST, sample sheet, parameters + sessionInfo |

## Minimum reproducibility information
* **Genome / annotation**: GRCm38 + GENCODE vM25 for every step (`config/reference_genome.txt`).
* **Software versions**: `config/software_versions.txt`.
* **Parameters and thresholds**: `config/parameters.txt`.
* **Transcript → gene → symbol**: tx2gene from the vM25 GTF; symbols from the GTF `gene_name`; duplicate symbols kept as
  separate genes; missing symbols fall back to the gene ID; Entrez (KEGG only) via org.Mm.eg.db (details in
  `config/reference_genome.txt`).
* **Custom signatures**: `metadata/gene_sets.csv`; score = mean per-gene z-score of VST expression.
* **Notes**: ND_10R_2 and ND_11L_2 were sequenced in a separate delivery (see `metadata/sample_metadata.csv`).

## Requirements
SLURM cluster; conda environment with STAR, samtools, salmon, cutadapt; TEtranscripts/TElocal; R ≥ 4.3 with the
packages listed in `config/software_versions.txt`. Raw FASTQ files are not included (deposited at GEO: [accession]).
