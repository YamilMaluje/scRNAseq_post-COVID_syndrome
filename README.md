# scRNAseq_post-COVID_syndrome

# 🧬 scRNA-seq Reveals Persistent Aberrant Differentiation of Nasal Epithelium Driven by TNFα and TGFβ in Post-COVID Syndrome

This repository contains the analysis code and workflows for our study:

**"scRNA-seq reveals persistent aberrant differentiation of nasal epithelium driven by TNFα and TGFβ in post-COVID syndrome"**  
📄 [bioRxiv DOI: 10.1101/2024.01.10.574801](https://doi.org/10.1101/2024.01.10.574801)

## 🧠 Overview

Our study reveals a pathway for persistent severe post-COVID syndrome (PCS) involving immune-mediated nasal tissue damage. We show that the TNFα-TGFβ axis drives aberrant differentiation in the nasal epithelium. Targeting this axis may restore epithelial homeostasis and offer a therapeutic strategy for PCS, akin to treatments in other chronic inflammatory conditions.

## 🔗 Data Availability

- **GEO Accession**: [GSE299529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299529)
- No example or downsampled data is included in this repository.

## 💻 Software & Tools

- Language: **R**
- Key package: [Seurat](https://satijalab.org/seurat/)
- Additional tools: Bash scripts for FASTQ QC and preprocessing

## ⚙️ Pipeline Overview

### Step 1: Quality Control of FASTQ Reads
bash
bash scripts/shell/QCFastq.sh
Tool: FASTQC

### Step 2: Read Trimming & Alignment
Adapter trimming with fastp

Alignment with kallisto (pseudoalignment)

### Step 3: Preprocessing & Filtering
Import Kallisto quantifications

Filter low-quality cells, normalize data

### Step 4: Downstream Analysis
Main downstream analysis (scRNA_Downstream_analysis_v2.1.R) 

Cell Annotation using canonical markers (Canonical_markers_Seurat_Annotation.R and Canonical_Markers_plots.R)

Differential Abundance Analysis (DA_Analysis_Cluster_v1.2.sh)

Differential Gene Expression (DGE_scRNAseq_Condition.R)

GSEA (Gene Set Enrichment Analysis) (GSEA.R)

Difference between the proprotion of cells in clusters (cell_proportions_v1.0.R)

Scripts are available under the scripts/ directory.

