#!/bin/bash

# ------------------------------------------------------------------
# Generate species-specific opsin gene lists.
#
# Hydra 105 and Hydra viridissima use protein IDs in the original
# opsin lists. These protein IDs are matched to their corresponding
# LOC IDs using the genome annotation files.
#
# Hydra vulgaris AEP and Hydra oligactis already use gene IDs.
# get_locs_simple.sh is used to generate the gene lists needed for
# the downstream analyses.
# ------------------------------------------------------------------

############################
# Hydra 105
############################

echo "Generating Hydra 105 opsin gene list..."

while read -r prot; do
    grep -m1 "$prot" GCF_022113875.1_Hydra_105_v3_genomic.gff \
    | grep -oP 'gene=LOC\d+' \
    | sed 's/gene=//'
done < hydra_105_opsin.txt > hydra_105_opsin_LOCs.txt

############################
# Hydra viridissima
############################

echo "Generating Hydra viridissima opsin gene list..."

while read -r prot; do
    grep -m1 "$prot" GCF_964215635.1_jhHydViri1.1_genomic.gff \
    | grep -oP 'gene=LOC\d+' \
    | sed 's/gene=//'
done < hydra_Viri_opsin.txt > hydra_Viri_LOCs.txt

############################
# Hydra vulgaris AEP
############################

echo "Generating Hydra vulgaris AEP opsin gene list..."

chmod +x get_locs_simple.sh

./get_locs_simple.sh \
    HVAEP_opsins.txt \
    HVAEP.GeneModels.gff3 \
    AEP_opsin_LOCs_derived.txt

############################
# Hydra oligactis
############################

echo "Generating Hydra oligactis opsin gene list..."

./get_locs_simple.sh \
    hydra_oli_opsins.txt \
    HOLI.GeneModels.gff3 \
    hydra_oli_opsin_LOCs.txt

echo "Finished generating opsin gene lists."
