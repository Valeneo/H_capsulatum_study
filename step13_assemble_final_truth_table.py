#!/usr/bin/env python3
"""
Step 13: Assemble Master Truth Table for the Histoplasma Diagnostic Kit
Tools: Python / Data Assembler

Biology & Logic:
- Merges the biophysical primer specifications (Forward/Reverse primers, Tm, amplicon sizes)
  with the multi-strain allele repeat matrix.
- Generates the definitive master laboratory blueprint:
  'TRUE_final_diagnostic_kit.tsv'
"""

import os

MATRIX_FILE = "results/final_primer_panel/top_20_strain_repeat_matrix.tsv"
PRIMER_MASTER = "results/final_primer_panel/universal_str_primers.tsv"
TRUE_FINAL_OUTPUT = "results/final_primer_panel/TRUE_final_diagnostic_kit.tsv"


def main():
    if not os.path.exists(MATRIX_FILE) or not os.path.exists(PRIMER_MASTER):
        print("❌ Error: Missing required input sheets. Please run Step 6 and Step 12 first.")
        return

    # 1. Load repeat matrix entries
    allele_matrix = {}
    with open(MATRIX_FILE, "r") as f:
        for line in f.readlines()[1:]:
            parts = line.strip().split("\t")
            if parts:
                locus_id = parts[0]
                motif = parts[1]
                alleles = parts[2:-1]  # all strain columns
                allele_matrix[locus_id] = {
                    "motif": motif,
                    "alleles_str": "|".join(alleles)
                }

    # 2. Match with real primer sequences
    print("📋 Compiling master diagnostic kit datasheet for Histoplasma capsulatum...")

    compiled_count = 0
    with open(TRUE_FINAL_OUTPUT, "w") as out:
        out.write("Locus_ID\tMotif\tForward_Primer_5to3\tReverse_Primer_5to3\tFwd_Tm\tRev_Tm\tExpected_Product_Size\tRepeat_Alleles_Per_Strain\n")

        with open(PRIMER_MASTER, "r") as pm:
            for line in pm.readlines()[1:]:
                parts = line.strip().split("\t")
                if not parts:
                    continue

                locus_id = parts[0]
                if locus_id in allele_matrix:
                    am = allele_matrix[locus_id]
                    fwd_seq = parts[4]
                    fwd_tm = parts[5]
                    rev_seq = parts[6]
                    rev_tm = parts[7]
                    p3_size = parts[8]

                    out.write(
                        f"{locus_id}\t{am['motif']}\t{fwd_seq}\t{rev_seq}\t{fwd_tm}\t{rev_tm}\t{p3_size}\t{am['alleles_str']}\n"
                    )
                    compiled_count += 1

    print(f"\n✨ Stage 13 Complete: Successfully assembled {compiled_count} validated diagnostic markers into:")
    print(f"📊 {TRUE_FINAL_OUTPUT}")


if __name__ == "__main__":
    main()
