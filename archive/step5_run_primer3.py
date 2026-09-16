#!/usr/bin/env python3
"""
Step 5: Automated Primer3 Input Formatting & High-Stringency Diagnostic Primer Design
Tools: primer3_core binary / BioPython

Biology & Logic:
- Loads the 100% cross-conserved 300bp-flanked STR sequences from Step 4.
- Formats Primer3 input blocks targeting the central STR repeat region (target at 301 bp).
- Enforces strict diagnostic PCR parameters:
  - Optimal Tm: 60.0°C (Range: 58.0°C - 62.0°C)
  - Max Tm difference between pairs: <= 2.0°C
  - Optimal GC content: 50.0% (Range: 40.0% - 60.0%)
  - Primer size: 18 - 25 bp (Optimal: 20 bp)
  - Expected Amplicon Product Size Range: 150 - 450 bp
  - Low self-complementarity (Max Self Any <= 8.0, Max Self End <= 3.0)
- Executes primer3_core and saves output to 'primer3/outputs/designed_primers.p3out'.
"""

import os
import subprocess
from Bio import SeqIO

CONSERVED_FASTA = "results/step3_conserved_sequences/conserved_all_genomes_flanks.fasta"
P3_DIR = "primer3/outputs"
os.makedirs(P3_DIR, exist_ok=True)
P3_INPUT_FILE = os.path.join(P3_DIR, "input_records.p3in")
P3_OUTPUT_FILE = os.path.join(P3_DIR, "designed_primers.p3out")


def main():
    if not os.path.exists(CONSERVED_FASTA):
        print(f"❌ Error: Conserved FASTA '{CONSERVED_FASTA}' not found. Please run Step 4 first.")
        return

    print("==========================================================================")
    print("🧬 Histoplasma capsulatum Universal Primer Design - Stage 5 (Primer3)")
    print("==========================================================================\n")

    records = list(SeqIO.parse(CONSERVED_FASTA, "fasta"))
    print(f"📋 Loaded {len(records)} conserved STR sequence blocks.")

    # 1. Format Primer3 input records
    formatted_count = 0
    with open(P3_INPUT_FILE, "w") as out:
        for rec in records:
            seq_str = str(rec.seq).upper()
            header_parts = rec.id.split("|")
            
            # Flank is exactly 300 bp left, target is the STR core
            try:
                coords = header_parts[2].split("-")
                repeat_len = int(coords[1]) - int(coords[0]) + 1
            except Exception:
                repeat_len = 30

            # Primer3 requires 0-based start for TARGET (300 = index of 301st base)
            out.write(f"SEQUENCE_ID={rec.id}\n")
            out.write(f"SEQUENCE_TEMPLATE={seq_str}\n")
            out.write(f"SEQUENCE_TARGET=300,{repeat_len}\n")
            out.write("PRIMER_TASK=generic\n")
            out.write("PRIMER_PICK_LEFT_PRIMER=1\n")
            out.write("PRIMER_PICK_INTERNAL_OLIGO=0\n")
            out.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
            out.write("PRIMER_OPT_SIZE=20\n")
            out.write("PRIMER_MIN_SIZE=18\n")
            out.write("PRIMER_MAX_SIZE=25\n")
            out.write("PRIMER_OPT_TM=60.0\n")
            out.write("PRIMER_MIN_TM=58.0\n")
            out.write("PRIMER_MAX_TM=62.0\n")
            out.write("PRIMER_PAIR_MAX_DIFF_TM=2.0\n")
            out.write("PRIMER_MIN_GC=40.0\n")
            out.write("PRIMER_OPT_GC_PERCENT=50.0\n")
            out.write("PRIMER_MAX_GC=60.0\n")
            out.write("PRIMER_MAX_SELF_ANY=8.00\n")
            out.write("PRIMER_MAX_SELF_END=3.00\n")
            out.write("PRIMER_PRODUCT_SIZE_RANGE=150-230 231-310 311-450\n")
            out.write("PRIMER_NUM_RETURN=5\n")
            out.write("=\n")
            formatted_count += 1

    print(f"📦 Formatted {formatted_count} locus blocks into Primer3 format: {P3_INPUT_FILE}")

    # 2. Run primer3_core
    print("🚀 Executing primer3_core calculation engine...")
    try:
        with open(P3_INPUT_FILE, "r") as in_f, open(P3_OUTPUT_FILE, "w") as out_f:
            subprocess.run(["primer3_core"], stdin=in_f, stdout=out_f, check=True)
        print(f"✅ Primer3 execution succeeded!")
        print(f"💾 Raw designs saved to: {P3_OUTPUT_FILE}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running primer3_core: {e}")
        return

    print("\n==========================================================================")
    print("🎉 STAGE 5 COMPLETE! You are ready for Step 6 (step6_parse_primer3_results.py).")
    print("==========================================================================")


if __name__ == "__main__":
    main()

