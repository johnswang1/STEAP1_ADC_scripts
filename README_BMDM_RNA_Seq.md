# RNA-seq Analysis: Dual-Species Transcriptomics of ADC Payloads

This repository contains the R analysis code accompanying the manuscript submitted to the Journal of Translational Medicine. The study characterizes the transcriptional response to antibody-drug conjugate (ADC) payloads (Exatecan and MMAE) in an in vitro co-culture model using dual-species (22Rv1.hSTEAP1 CM cells + mouse BMDM) bulk RNA-seq.

---

## Experimental Design

| Group | n | Description |
|---|---|---|
| Isotype | 3 (A1–A3) | Isotype control |
| WT | 3 (A4–A6) | Untreated co-culture |
| WT_Exatecan | 3 (A7–A9) | ADC with Exatecan payload |
| WT_MMAE | 3 (A10–A12) | ADC with MMAE payload |

Samples contain both human and mouse cell transcripts from in vitro co-culture. Reads were aligned to a combined human (GRCh38) + mouse (GRCm39) reference to deconvolve species of origin. All analyses were performed on mouse BMDM reads only.

---

## Repository Structure

```
.
├── metadata.csv                  # Sample-to-condition mapping
├── RNA_Seq code.Rmd              # Code file 1 (primary DESeq2 pipeline)
├── RNA_Seq code 2.Rmd            # Code file 2 (heatmaps, GO, Reactome, volcanos)
├── reference/                    # Reference genomes and transcriptomes
│   ├── Homo_sapiens.GRCh38.*    # Human GRCh38 (Ensembl release 115)
│   ├── Mus_musculus.GRCm39.*    # Mouse GRCm39 (Ensembl release 115)
│   ├── human_mouse_combined_txome.gtf   # Combined GTF used for tximport
│   ├── human_mouse_transcriptome.fa     # Combined FASTA for Salmon index
│   ├── salmon_human_mouse_index/        # Pre-built Salmon index
│   └── decoys.txt               # Decoy sequences for selective alignment
├── salmon/                       # Salmon quantification output
│   └── {SampleID}/quant.sf      # Per-sample transcript quantification
├── results/                      # Output tables and figures
└── qc/                           # FastQC / MultiQC reports
```

---

## Upstream Processing (Terminal)

The following steps were run in the terminal prior to R analysis. Salmon version 1.10 was used.

### 1. Build combined human + mouse Salmon index

Species prefixes (`human_`, `mouse_`) were added to all transcript and chromosome names to keep species identifiable after joint quantification.

```bash
salmon index \
  --transcripts human_mouse_transcriptome.fa \
  --decoys decoys.txt \
  --index reference/salmon_human_mouse_index \
  --gencode \
  -p 8
```

### 2. Quantify each sample with Salmon

```bash
for SAMPLE in A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12; do
  salmon quant \
    --index reference/salmon_human_mouse_index \
    --libType A \
    --validateMappings \
    --gcBias \
    --numBootstraps 30 \
    -r fastq/${SAMPLE}.fastq.gz \
    --output salmon/${SAMPLE}
done
```

> **Note:** `--libType A` lets Salmon detect library strandedness automatically. Adjust `-r` to `-1`/`-2` flags for paired-end data.

---

## R Analysis

All downstream analysis is performed in R (version ≥ 4.3) using Bioconductor packages.

### Installation

Open `RNA_Seq code 2.Rmd` in RStudio and run the `setup-install` chunk once (set `eval=TRUE` temporarily):

```r
BiocManager::install(c(
  "DESeq2", "apeglm", "tximport", "rtracklayer",
  "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi",
  "EnhancedVolcano", "pheatmap", "ReactomePA", "enrichplot"
))
install.packages(c("tidyverse", "ggrepel", "stringr", "forcats"))
```

### Running the analysis

1. Clone this repository and open the project root in RStudio.
2. Confirm `project_dir <- "."` in the `set-paths` chunk points to the repo root.
3. Knit `RNA_Seq code 2.Rmd` or run chunks sequentially.

### Analysis steps (`RNA_Seq code 2.Rmd`)

| Chunk | Description |
|---|---|
| `tx2gene-and-tximport` | Import Salmon counts via tximport; build tx-to-gene map from combined GTF |
| `split-human-mouse` | Separate human and mouse count matrices; compute species composition per sample |
| `deseq2-mouse` | DESeq2 differential expression on mouse counts |
| `mhc1-heatmap` | Heatmap of MHC-I antigen presentation genes |
| `inflammatory-heatmap` | Heatmap of inflammatory macrophage-associated genes |
| `go-enrichment` | GO Biological Process enrichment (WT_MMAE vs. WT_Exatecan) |
| `reactome-exa` / `reactome-exa-tiff` | Reactome pathway enrichment — Exatecan-enriched pathways |
| `supp-reactome-mmae` | Supplemental: Reactome — MMAE-enriched pathways |
| `supp-volcano-wt-vs-isotype-complement` | Supplemental: WT vs. Isotype volcano highlighting complement genes |

---

## Key Outputs

| File | Description |
|---|---|
| `results/mouse_WT_Exatecan_vs_WT_MMAE.csv` | DESeq2 results table |
| `results/species_fraction_summary.csv` | Human vs. mouse read fractions per sample |
| `results/MHC1_antigen_presentation_heatmap.tiff` | MHC-I antigen presentation gene heatmap |
| `results/Inflammatory_macrophage_genes_heatmap.tiff` | Inflammatory macrophage gene heatmap |
| `results/Reactome_WT_Exatecan_vs_MMAE_with_IDs.tiff` | Main Reactome figure |
| `results/Supp_Complement_dotplot_all_conditions.tiff` | Complement gene expression |

---

## Reference Genomes

| Species | Assembly | Ensembl Release | Source |
|---|---|---|---|
| Human | GRCh38 | 115 | Ensembl |
| Mouse | GRCm39 | 115 | Ensembl |

---

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon publication) and the following tools:

- **Salmon:** Patro et al. (2017) *Nature Methods* 14, 417–419
- **DESeq2:** Love et al. (2014) *Genome Biology* 15, 550
- **tximport:** Soneson et al. (2015) *F1000Research* 4, 1521
- **clusterProfiler:** Wu et al. (2021) *The Innovation* 2(3), 100141
- **ReactomePA:** Yu & He (2016) *Molecular BioSystems* 12, 477–479
