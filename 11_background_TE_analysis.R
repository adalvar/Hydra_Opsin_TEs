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

## ClassI, ClassII, and unknown stay fixed during this rotation;
## the rest shift down one position.

chain <- c(
  "DNA", "LINE", "LTR", "MITE", "nLTR",
  "nMITE", "PLE", "RC", "SINE"
)

stopifnot(all(chain %in% names(spectral_colors)))

shifted_colors <- spectral_colors

shifted_colors[chain] <- spectral_colors[
  c(
    "SINE", "DNA", "LINE", "LTR", "MITE",
    "nLTR", "nMITE", "PLE", "RC"
  )
]

## ---- Swap ClassII and nMITE colors ----

swap_tmp <- shifted_colors["ClassII"]

shifted_colors["ClassII"] <- shifted_colors["nMITE"]

shifted_colors["nMITE"] <- swap_tmp

############################################################
# Plot mean background TE coverage
############################################################

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
# Calculate background coverage variation
############################################################

# Number of gene x region units.
# Missing TE classes within observed units count as zero coverage.

n_units <- n_distinct(bg_df$gene, bg_df$region)

# Mean per class across all units (stack segments).

bg_avg <- bg_df %>%
  group_by(TE_class) %>%
  summarise(
    pct_coverage = sum(pct_coverage) / n_units,
    .groups = "drop"
  ) %>%
  mutate(TE_class = fct_reorder(TE_class, pct_coverage))

# Total TE coverage per gene x region (one number per unit).

unit_totals <- bg_df %>%
  group_by(gene, region) %>%
  summarise(total = sum(pct_coverage), .groups = "drop")

# Mean total coverage equals mean(unit_totals$total).

total_mean <- sum(bg_avg$pct_coverage)
total_sd <- sd(unit_totals$total)
total_ci <- qt(0.975, n_units - 1) * total_sd / sqrt(n_units)

# Mean +/- one standard deviation.

c(
  lower = total_mean - total_sd,
  upper = total_mean + total_sd
)

# Approximate 95% confidence interval edges.
# Assumes independent gene x region units.

c(
  lower = total_mean - total_ci,
  upper = total_mean + total_ci
)

# Central 80% interval of background coverage.
# The 10th percentile is the coverage below which 10% of units fall.
# The 90th percentile is the coverage above which 10% of units fall.

quantile(unit_totals$total, c(0.10, 0.90))

############################################################
# Expected output
############################################################

# Mean background TE coverage plot.
#
# Mean total background TE coverage for comparison with
# opsin-associated TE coverage.
#
# Mean +/- one standard deviation.
#
# Approximate 95% confidence interval for the mean.
#
# Central 80% interval of background coverage.

message("Finished calculating background TE coverage.")
