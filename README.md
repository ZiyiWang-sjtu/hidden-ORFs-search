# Hidden ORFs Search

A computational pipeline for identifying potentially expressible hidden open reading frames (ORFs) in plasmid sequences.

## Overview

DNA sequences can encode proteins in six possible reading frames, which means plasmids may contain unintended open reading frames outside the intended coding sequence. Some of these hidden ORFs may fall within transcriptional units and potentially be expressed.

This project provides a computational pipeline to identify hidden ORFs in plasmid sequences and evaluate their potential expression based on regulatory elements and their positional relationship to the target gene.

## Requirements

The pipeline was tested with the following environment:

- Python 3.9
- Biopython 1.85
- pandas 2.3
- numpy 2.0
- BLAST+ 2.17
- tqdm 4.67

### R libraries required

The following R packages are required for data analysis and visualization:

- tidyverse
- scales
- readr

## 1. Filter Plasmids Containing Target Genes

### Input files

- JSON file containing plasmid information, which is provided by the addgene team.
- FASTA file containing the coding sequences (CDS) of target genes of interest.

### Example command

```bash
# Filter plasmids containing target genes
python3 plasmid_json_filter.py --input addgene-plasmids-sequences.json --fasta target_gene.fa --output filtered_target_plasmids.json 
```
## 2. Extract Features from SnapGene

This script extracts specific genetic features (e.g., promoters, polyA signals, terminators, LTR) from SnapGene files. It standardizes feature types and compiles them into a CSV library for downstream analysis.

### Input files

Folders containing .dna files exported from SnapGene with annotated features. Place files into folders based on feature type.

### Example command

```bash
# Extract features from Snapgene
python3 SnapGene_Feature_Extractor.py --input features_terminator_from_snapgene features_promoter_from_snapgene features_polya_from_snapgene 
```

We additionally manually curated the transcriptional regulatory elements that are functional in mammalian cells. The curated list is provided in mammalian_feature_library.csv.

## 3. Hidden ORFs Search

This script performs a comprehensive analysis of plasmid sequences, mapping regulatory features and scanning for expressible ORFs.

### 3.1 Input and example command

- --json: plasmid JSON file from plasmid_json_filter.py
- --fa: target gene FASTA file
- --features: feature library CSV from SnapGene_Feature_Extractor.py

Optional parameters:

- --identity: minimum sequence identity (%) for BLAST hits (default 96%)

- --min_aa: minimum amino acid length for ORFs (default 100)

- --max_prom: maximum allowed distance from upstream promoter (default 10,000 bp)

- --max_poly: maximum allowed distance to downstream polyA (default 10,000 bp)

- --out: output CSV file name

- --workers: number of parallel processes to use


```bash
# Hidden ORFs Search
python3 hidden_orf_search.py --json filtered_target_plasmids.json --fa target_gene.fa --features feature_library.csv --min_aa 100 --max_prom 10000 --max_poly 10000 --out results.csv
```

### 3.2 Preprocessing data and target gene mapping 

### 3.3 ORF scanning
ORFs are scanned across all frames on both strands of the circular plasmid, and only ORFs longer than min_aa are retained.

### 3.4 ORF classification relative to target gene
Each ORF is classified based on its overlap with the target gene:
- inside_target_gene: Fully contained within the GOI.
- spanning_target_gene: Overlaps the entire GOI.
- partial_target_gene: Partial overlap.
- outside_target_gene: No overlap.

> [!NOTE]
> **This pipeline supports two analysis modes:**
>
> **1. Hidden ORFs (DNA-level analysis)**  
> This mode identifies all ORFs that exist at the DNA sequence level based solely on canonical translation signals (an **ATG** start codon and an in-frame stop codon). No transcriptional regulatory elements are considered in this analysis.
>
> **2. Putatively Expressible Hidden ORFs (Transcription-aware analysis)**  
> This mode further evaluates whether a Hidden ORF is potentially transcribable by requiring both an upstream **transcription initiation element** and a downstream **transcription termination element** on the same strand and in the correct transcriptional orientation.

### 3.5 Feature mapping and distance calculation

### 3.6 Expressibility and risk assessment
