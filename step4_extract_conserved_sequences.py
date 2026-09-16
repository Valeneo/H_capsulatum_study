#!/usr/bin/env python3
"""
Step 4: Extract Conserved Sequence Blocks for Primer3 Input
Tools: BioPython (SeqIO)

Biology & Logic:
- Loads the verified universally conserved STR locus IDs identified in Step 3.
- Extracts the complete 300bp-flanked sequence records from the query consensus FASTA.
- Formats and writes the FASTA sequence blocks ready for high-stringency Primer3 primer design.
"""

import os
import glob
from Bio import SeqIO

CONSERVED_META = "results/step2_conserved_loci/conserved_loci_metadata.txt"
FLANK_DIR = "flanking_sequences"
OUTPUT_DIR = "results/step3_conserved_sequences"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FASTA = os.path.join(OUTPUT_DIR, "conserved_all_genomes_flanks.fasta")


def main():
    if not os.path.exists(CONSERVED_META):
        print(f"❌ Error: Conserved metadata '{CONSERVED_META}' not found. Please run Step 3 first.")
        return

    # 1. Load conserved locus IDs
    conserved_ids = set()
    with open(CONSERVED_META, "r") as f:
        for line in f.readlines()[1:]:  # Skip header
            parts = line.strip().split("\t")
            if parts:
                conserved_ids.add(parts[0])

    print(f"📋 Loaded {len(conserved_ids)} target conserved STR locus IDs from metadata.")

    # 2. Find available query flank FASTA files
    flank_files = glob.glob(os.path.join(FLANK_DIR, "*_consensus_flanks.fasta"))
    if not flank_files:
        print(f"❌ Error: No flank FASTA files found in '{FLANK_DIR}/'. Please run Step 2 first.")
        return

    # Use the primary query flank file (or first available)
    query_fasta = flank_files[0]
    print(f"🧬 Extracting sequences from primary reference: '{os.path.basename(query_fasta)}'...")

    extracted_count = 0
    with open(OUTPUT_FASTA, "w") as out_fa:
        for record in SeqIO.parse(query_fasta, "fasta"):
            if record.id in conserved_ids:
                out_fa.write(f">{record.id}\n{record.seq}\n")
                extracted_count += 1

    print(f"\n✨ Stage 4 Complete: Successfully isolated {extracted_count} universally conserved sequence blocks.")
    print(f"📁 Primer3-ready input file saved at:\n   {OUTPUT_FASTA}")


if __name__ == "__main__":
    main()
