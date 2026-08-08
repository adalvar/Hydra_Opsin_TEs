# ------------------------------------------------------------------
# Calculate mean transposable element coverage across background genes.
#
# In this step, transposable element coverage surrounding background
# protein-coding genes is summarized by TE superfamily. The resulting
# average coverage values are used as the genomic background for
# comparison with opsin-associated TE coverage.
# ------------------------------------------------------------------

############################################################
# Load required packages
############################################################

library(tidyverse)
library(RColorBrewer)

############################################################
# Set working directory
############################################################

# Set the working directory to the folder containing the
# background TE coverage tables.

setwd("<PATH_TO_PROJECT_DIRECTORY>")

############################################################
# Input file
############################################################

# Replace <SPECIES> with the genome being analyzed
# (e.g., AEP, Oli, 105, Viri).

input_file <- "<SPECIES>_background_TE_coverage_by_class.tsv"

############################################################
# R analysis script
############################################################

# Paste the complete R analysis script below.

bg_df <- read_tsv(input_file)

bg_avg <- bg_df %>%
  group_by(TE_class) %>%
  summarise(pct_coverage = mean(pct_coverage), .groups = "drop")

## ---- Spectral palette, interpolated to cover all TE_class levels ----

te_levels <- sort(unique(bg_df$TE_class))

spectral_colors <- setNames(
  colorRampPalette(brewer.pal(11, "Spectral"))(length(te_levels)),
  te_levels
)

## ---- Rotate colors through TE superfamily categories ----

## ClassI, ClassII, and unknown stay fixed; the rest shift down one position

chain <- c("DNA", "LINE", "LTR", "MITE", "nLTR", "nMITE", "PLE", "RC", "SINE")

stopifnot(all(chain %in% names(spectral_colors)))

shifted_colors <- spectral_colors

shifted_colors[chain] <- spectral_colors[
  c("SINE", "DNA", "LINE", "LTR", "MITE",
    "nLTR", "nMITE", "PLE", "RC")
]

## ---- Swap ClassII and nMITE colors ----

swap_tmp <- shifted_colors["ClassII"]

shifted_colors["ClassII"] <- shifted_colors["nMITE"]

shifted_colors["nMITE"] <- swap_tmp

ggplot(bg_avg,
       aes(x = "Background",
           y = pct_coverage,
           fill = TE_class)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = shifted_colors) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Mean % TE coverage",
    fill = "TE class"
  ) +
  theme_minimal()

sum(bg_avg$pct_coverage)

############################################################
# Expected output
############################################################

#
# Mean background TE coverage plot
#
# The mean background TE coverage value is used as the
# reference line in the opsin TE coverage figures.
#

message("Finished calculating background TE coverage.")
