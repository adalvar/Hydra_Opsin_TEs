# _Hydra_ Opsin TE Analysis

Code and analysis pipelines for investigating transposable element composition surrounding *Hydra* opsin genes.


# Overview

This repository contains the scripts and workflow used to characterize transposable element (TE) composition surrounding opsin genes in multiple *Hydra* genomes.


# Manuscript

This repository accompanies an ongoing study investigating transposable element composition surrounding opsin genes across multiple *Hydra* species.

*Citation will be added upon publication.*


# Platforms Used

- UCSC Hummingbird Computing Cluster (command line)
- RStudio


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


# Methodology

1. Download genome assemblies and corresponding genome annotation files.

2. Annotate transposable elements using RepeatMasker with a custom Medusozoa repeat library.

3. Generate species-specific opsin gene lists from genome annotation files.

4. Create sorted BED files containing genomic coordinates for all opsin genes.

5. Index each genome using `samtools faidx` and generate genome size files for BEDTools.

6. Convert RepeatMasker output (`.out`) into sorted BED files containing TE coordinates and classifications.

7. Generate 50 kb upstream and downstream flanking regions surrounding each opsin gene using BEDTools.

8. Calculate total TE coverage within each flanking region.

9. Separate TEs by superfamily, resolve overlapping annotations, and calculate per-class TE coverage.

10. Export TE coverage summaries as tab-delimited (`.tsv`) files.

11. Generate visualizations of TE composition using custom R scripts.

12. Calculate genome-wide background TE coverage for comparison with opsin-associated regions.


# Repository Contents

```
scripts/
    Shell, Python, and SLURM scripts used throughout the analysis pipeline

R/
    R scripts used for visualization and downstream analyses

example_files/
    Example input files and templates

docs/
    Additional documentation
```

# Output Files

The pipeline produces:

- TE coverage summaries (`.tsv`)
- Total upstream TE coverage
- Total downstream TE coverage
- Publication-quality PDF figures
- Background TE coverage summaries


# Data Availability

Genome assemblies and genome annotation files were obtained from NCBI GenBank. Repeat annotations were generated using RepeatMasker with a custom Medusozoa repeat library.


# Acknowledgements

This pipeline was developed in the Macias-Muñoz Lab at the University of California, Santa Cruz.

