#!/usr/bin/env python3

"""
Create sorted BED files containing opsin gene coordinates.

This script matches species-specific opsin IDs to their corresponding
genomic coordinates in a GFF3 annotation file and generates a sorted
BED file for downstream TE analyses.

Usage:
    python3 make_opsin_bed.py <gff3_file> <opsin_id_list> <output.bed>
"""

import re
import sys
from pathlib import Path

def parse_gff3_genes(gff3_path):
    """Return dict: lowercased gene ID -> (scaffold, start, end, strand)."""
    gene_lookup = {}
    id_pattern = re.compile(r"ID=([^;]+)")

    with open(gff3_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != "gene":
                continue

            scaffold = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            match = id_pattern.search(attributes)
            if not match:
                continue

            gene_id = match.group(1)
            gene_lookup[gene_id.lower()] = (scaffold, start, end, strand)

    return gene_lookup


def load_opsin_ids(opsin_list_path):
    ids = []
    with open(opsin_list_path) as fh:
        for line in fh:
            gene_id = line.strip()
            if gene_id:
                ids.append(gene_id)
    return ids


def build_bed_rows(gene_lookup, opsin_ids):
    rows = []
    missing = []

    for opsin_id in opsin_ids:
        entry = gene_lookup.get(opsin_id.lower())
        if entry is None:
            missing.append(opsin_id)
            continue

        scaffold, start, end, strand = entry
        bed_start = start - 1  # GFF3 is 1-based; BED start is 0-based
        rows.append((scaffold, bed_start, end, opsin_id, ".", strand))

    return rows, missing


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: python3 make_opsin_bed.py <gff3_file> <opsin_id_list> <output.bed>",
            file=sys.stderr,
        )
        sys.exit(1)

    gff3_path = Path(sys.argv[1])
    opsin_list_path = Path(sys.argv[2])
    output_path = Path(sys.argv[3])

    gene_lookup = parse_gff3_genes(gff3_path)
    opsin_ids = load_opsin_ids(opsin_list_path)
    rows, missing = build_bed_rows(gene_lookup, opsin_ids)

    rows.sort(key=lambda r: (r[0], r[1]))

    with open(output_path, "w") as out:
        for scaffold, start, end, gene_id, score, strand in rows:
            out.write(f"{scaffold}\t{start}\t{end}\t{gene_id}\t{score}\t{strand}\n")

    print(f"Wrote {len(rows)} opsin entries to {output_path}")

    if missing:
        print(f"\nWARNING: {len(missing)} opsin ID(s) not found in the GFF3 file:")
        for gene_id in missing:
            print(f"  {gene_id}")


if __name__ == "__main__":
    main()
