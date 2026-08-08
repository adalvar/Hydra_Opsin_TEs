# ------------------------------------------------------------------
# Analyze transposable element coverage surrounding opsin genes.
#
# In this step, annotated TE coverage tables are imported into R,
# summarized by TE superfamily, and visualized using stacked bar
# plots. Separate plots are generated for upstream, downstream,
# and combined flanking regions surrounding opsin genes.
#
# Opsin genes are displayed according to the phylogenetic clade
# order defined by Rachel P.'s opsin phylogenies to maintain a
# consistent biological ordering across all figures.
# ------------------------------------------------------------------

############################################################
# Load required packages
############################################################

library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)

############################################################
# Set working directory
############################################################

# Set the working directory to the folder containing the
# annotated TE coverage tables.

setwd("<PATH_TO_PROJECT_DIRECTORY>")

############################################################
# Input file
############################################################

# Replace <SPECIES> with the genome being analyzed
# (e.g., AEP, Oli, 105, Viri).

input_file <- "<SPECIES>_TE_coverage_by_class_out.tsv"

############################################################
# R analysis script
############################################################

##Hydra viridissima (Hviri)
df <- read_tsv("Viri_TE_coverage_by_class_out.tsv", show_col_types = FALSE) %>%
  mutate(
    region = factor(region, levels = c("upstream", "downstream")),
    clade = str_extract(Gene_ID_phylo, "(?<=Ops_)[A-Za-z]+")
  ) %>%
  arrange(clade, gene) %>%
  mutate(gene = factor(gene, levels = unique(gene)))

## ---- Collapse upstream + downstream into total coverage per gene/class ----
df_total <- df %>%
  group_by(gene, Gene_ID_phylo, TE_class) %>%
  summarise(pct_coverage = mean(pct_coverage), .groups = "drop") %>%
  mutate(region = "total")

## ---- Spectral palette, interpolated to cover all TE_class levels ----
te_levels <- sort(unique(df$TE_class))
spectral_colors <- setNames(
  colorRampPalette(brewer.pal(11, "Spectral"))(length(te_levels)),
  te_levels
)

## ---- Rotate colors through TE superfamily categories ----
## ClassI, ClassII, and unknown stay fixed; the rest shift down one position
chain <- c("DNA", "LINE", "LTR", "MITE", "nLTR", "nMITE", "PLE", "RC", "SINE")
stopifnot(all(chain %in% names(spectral_colors)))

shifted_colors <- spectral_colors
shifted_colors[chain] <- spectral_colors[c("SINE", "DNA", "LINE", "LTR", "MITE", "nLTR", "nMITE", "PLE", "RC")]

## ---- Swap ClassII and nMITE colors ----
swap_tmp <- shifted_colors["ClassII"]
shifted_colors["ClassII"] <- shifted_colors["nMITE"]
shifted_colors["nMITE"] <- swap_tmp

## ---- Single plot: genes on x (labeled by full clade ID), stacked % TE coverage on y, faceted by region ----
p <- ggplot(df, aes(x = Gene_ID_phylo, y = pct_coverage, fill = TE_class)) +
  geom_col() +
  facet_wrap(~region, ncol = 1) +
  scale_x_discrete(limits = c("Hviri_Ops_CII_06a","Hviri_Ops_CII_06b",
                              "Hviri_Ops_CIII_01","Hviri_Ops_CIII_02a","Hviri_Ops_CIII_02b","Hviri_Ops_CIII_03",
                              "Hviri_Ops_CIII_04a","Hviri_Ops_CIII_04b","Hviri_Ops_CIII_05","Hviri_Ops_CIII_06",
                              "Hviri_Ops_CIV_01","Hviri_Ops_CIV_02","Hviri_Ops_CIV_03","Hviri_Ops_CIV_04b",
                              "Hviri_Ops_CIV_05a","Hviri_Ops_CIV_05b","Hviri_Ops_CIV_06","Hviri_Ops_CIV_07","Hviri_Ops_CIV_08",
                              "Hviri_Ops_CIV_09a","Hviri_Ops_CIV_09b","Hviri_Ops_CIV_10")) +
  scale_fill_manual(values = shifted_colors) +
  ## TODO: replace with Hviri's mean background TE coverage (sum(bg_avg$pct_coverage) from the Viri background script)
  geom_hline(yintercept = 45.18337, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Opsin gene (Clade ID)", y = "% TE coverage",
    fill = "TE superfamily",
    title = "TE superfamily coverage in 50kb flanks of opsin genes (Hviri)"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "right")

## ---- Total-only plot ----
p_total <- ggplot(df_total, aes(x = Gene_ID_phylo, y = pct_coverage, fill = TE_class)) +
  geom_col() +
  scale_x_discrete(limits = c("Hviri_Ops_CII_06a","Hviri_Ops_CII_06b",
                              "Hviri_Ops_CIII_01","Hviri_Ops_CIII_02a","Hviri_Ops_CIII_02b","Hviri_Ops_CIII_03",
                              "Hviri_Ops_CIII_04a","Hviri_Ops_CIII_04b","Hviri_Ops_CIII_05","Hviri_Ops_CIII_06",
                              "Hviri_Ops_CIV_01","Hviri_Ops_CIV_02","Hviri_Ops_CIV_03","Hviri_Ops_CIV_04b",
                              "Hviri_Ops_CIV_05a","Hviri_Ops_CIV_05b","Hviri_Ops_CIV_06","Hviri_Ops_CIV_07","Hviri_Ops_CIV_08",
                              "Hviri_Ops_CIV_09a","Hviri_Ops_CIV_09b","Hviri_Ops_CIV_10")) +
  scale_fill_manual(values = shifted_colors) +
  ## TODO: replace with Hviri's mean background TE coverage (sum(bg_avg$pct_coverage) from the Viri background script)
  geom_hline(yintercept = 45.18337, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Opsin gene (Clade ID)", y = "% TE coverage (100kb flank total)",
    fill = "TE superfamily",
    title = "TE superfamily coverage in 100kb flanks (upstream + downstream) of opsin genes (Hviri)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "right"
  )

## ---- Save ----
ggsave("Viri_TE_coverage_by_Clade.pdf", p, width = 14, height = 8)
ggsave("Viri_TE_coverage_by_Clade_total.pdf", p_total, width = 14, height = 8)
p
p_total

############################################################
# Expected output
############################################################

#
# <SPECIES>_TE_coverage_by_Clade.pdf
# <SPECIES>_TE_coverage_by_Clade_total.pdf
#
# These figures summarize transposable element coverage
# surrounding opsin genes by TE superfamily and are used
# for downstream comparative analyses.
#

message("Finished generating opsin TE coverage figures.")
