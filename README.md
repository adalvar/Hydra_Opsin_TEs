# *Hydra* Opsin TE Analysis

Code and analysis pipelines used to characterize transposable element (TE) composition surrounding opsin genes across multiple *Hydra* genomes.

---

# Overview

This repository contains the scripts and workflow used to characterize transposable element composition surrounding opsin genes in multiple *Hydra* species. The pipeline includes repeat annotation, generation of genomic feature files, calculation of TE coverage surrounding opsin and background genes, mapping of phylogenetic gene names, and downstream visualization in R.

---

# Manuscript

This repository accompanies an ongoing study investigating transposable element composition surrounding opsin genes across multiple *Hydra* species.

*Citation will be added upon publication.*

---

# Platforms Used

- UCSC Hummingbird Computing Cluster (Linux command line)
- RStudio

---

# Packages and Dependencies

## Genome Annotation

- RepeatModeler v2.0.5
- RepeatMasker
- DeepTE
- BEDTools
- SAMtools

## Python

- Python 3

## R

- tidyverse
- ggplot2
- dplyr
- RColorBrewer

## Command Line Utilities

- awk
- grep
- sed
- sort

---

# Workflow

1. Generate species-specific opsin gene lists from genome annotation files.
2. Annotate transposable elements using RepeatMasker with a custom Medusozoa repeat library.
3. Create sorted BED files containing genomic coordinates for opsin genes.
4. Generate genome index and chromosome size files using `samtools faidx`.
5. Convert RepeatMasker output (`.out`) into sorted BED files containing TE coordinates and classifications.
6. Calculate transposable element coverage within 50 kb upstream and downstream flanking regions surrounding opsin genes.
7. Generate sorted BED files containing background protein-coding genes.
8. Calculate genome-wide background TE coverage using the same workflow.
9. Map species-specific opsin gene identifiers to shared phylogenetic gene names.
10. Generate publication-quality figures and summary statistics using custom R scripts.

---

# Repository Contents

```
01_generate_opsin_gene_lists.sh
02_repeatmasker.sh
03_make_opsin_bed.py
04_create_genome_info.sh
05_repeatmasker_to_bed.sh
06_opsin_TE_flanking.sh
07_make_background_gene_bed.sh
08_background_TE_flanking.sh
09_map_phylogenetic_gene_names.sh
10_opsin_TE_analysis.R
11_background_TE_analysis.R
README.md
```

---

# Input Files

The workflow requires:

- Genome assembly FASTA files
- Genome annotation files (`.gff` or `.gff3`)
- Species-specific opsin gene lists
- RepeatMasker output files (`.out`)
- `Opsin_Provenance_Master.csv`
- Custom Medusozoa repeat library

---

# Output Files

The pipeline generates:

- Sorted opsin BED files
- Genome index and chromosome size files
- Sorted RepeatMasker BED files
- Opsin TE coverage tables (`*_TE_coverage_by_class.tsv`)
- Background TE coverage tables (`*_background_TE_coverage_by_class.tsv`)
- Annotated TE coverage tables (`*_TE_coverage_by_class_out.tsv`)
- Publication-quality PDF figures

---

# Data Availability

Genome assemblies and genome annotation files were obtained from NCBI GenBank.

Repeat annotations were generated using RepeatMasker with a custom Medusozoa repeat library.

The `Opsin_Provenance_Master.csv` file, used to map species-specific opsin gene identifiers to shared phylogenetic gene names, was provided as part of the project resources.

---

# Acknowledgements

This computational workflow was developed in the Cnido Lab at the University of California, Santa Cruz.
