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
import sys
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
        sys.exit(1)

    # 1. Load conserved locus IDs
    conserved_ids = set()
    with open(CONSERVED_META, "r") as f:
        for line in f.readlines()[1:]:  # Skip header
            parts = line.strip().split("\t")
            if parts:
                conserved_ids.add(parts[0])

    if not conserved_ids:
        print(f"❌ Error: No conserved loci found in '{CONSERVED_META}'.")
        sys.exit(1)

    print(f"📋 Loaded {len(conserved_ids)} target conserved STR locus IDs from metadata.")

    # 2. Find matching query flank FASTA file
    ref_file = os.path.join(FLANK_DIR, "Hcaps_G217B_Reference_consensus_flanks.fasta")
    if os.path.exists(ref_file):
        query_fasta = ref_file
    else:
        # Fallback to any consensus flanks matching the prefix
        sample_id = next(iter(conserved_ids))
        prefix = sample_id.split("|")[0].replace("STR_HC_", "")
        matches = glob.glob(os.path.join(FLANK_DIR, f"*{prefix}*_consensus_flanks.fasta"))
        if matches:
            query_fasta = matches[0]
        else:
            all_flanks = glob.glob(os.path.join(FLANK_DIR, "*_consensus_flanks.fasta"))
            if not all_flanks:
                print(f"❌ Error: No flank FASTA files found in '{FLANK_DIR}/'.")
                sys.exit(1)
            query_fasta = all_flanks[0]

    print(f"🧬 Extracting sequences from primary reference: '{os.path.basename(query_fasta)}'...")

    extracted_count = 0
    with open(OUTPUT_FASTA, "w") as out_fa:
        for record in SeqIO.parse(query_fasta, "fasta"):
            if record.id in conserved_ids:
                out_fa.write(f">{record.id}\n{record.seq}\n")
                extracted_count += 1

    if extracted_count == 0:
        print(f"❌ Error: Extracted 0 matching sequences. Check ID formatting between {query_fasta} and {CONSERVED_META}.")
        sys.exit(1)

    print(f"\n✨ Stage 4 Complete: Successfully isolated {extracted_count} universally conserved sequence blocks.")
    print(f"📁 Primer3-ready input file saved at:\n   {OUTPUT_FASTA}")


if __name__ == "__main__":
    main()
