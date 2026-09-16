#!/usr/bin/env python3
"""
Step 7: Tri- and Tetra-Nucleotide Motif Filtering & Multiplex Size Binning
Tools: Python / Regular Expressions

Biology & Logic:
- Tri- (3bp) and Tetra- (4bp) nucleotide repeats produce significantly cleaner PCR amplicon bands
  without the polymerase stutter / slip artifacts typical of dinucleotides.
- Higher repeat copy numbers correlate with higher genetic instability and higher discriminatory power.
- Bins candidates into 3 non-overlapping amplicon size brackets:
  1. Small : 150 - 230 bp
  2. Medium: 231 - 310 bp
  3. Large : 311 - 450 bp
- Enables multiplex PCR where all 3 markers can be resolved simultaneously on a single gel/capillary channel.
"""

import os
import re

PRIMER_TABLE = "results/final_primer_panel/universal_str_primers.tsv"
OUTPUT_DIR = "results/final_primer_panel"
os.makedirs(OUTPUT_DIR, exist_ok=True)
TOP_MARKERS_FILE = os.path.join(OUTPUT_DIR, "top_universal_candidates.tsv")


def main():
    if not os.path.exists(PRIMER_TABLE):
        print(f"❌ Error: Primer table '{PRIMER_TABLE}' not found. Please run Step 6 first.")
        return

    candidates = []

    with open(PRIMER_TABLE, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            locus_id = parts[0]
            chrom = parts[1]
            coords = parts[2]
            motif = parts[3]
            fwd_seq = parts[4]
            fwd_tm = float(parts[5])
            rev_seq = parts[6]
            rev_tm = float(parts[7])
            amp_size = int(parts[8])

            # Filter strictly for Tri- and Tetra-nucleotide repeat units
            # Matches motifs like (GAA)12, (AGTC)8, or compound repeats like (CTG)8(TTG)14
            units = re.findall(r'\(([A-Za-z]{3,4})\)(\d+)', motif)

            if units:
                # Calculate total repeat unit count as proxy for mutation/polymorphism potential
                total_repeat_units = sum(int(count) for unit, count in units)
                # Ensure minimum 5 repeat units for adequate polymorphism
                if total_repeat_units >= 5:
                    candidates.append({
                        "locus_id": locus_id,
                        "chrom": chrom,
                        "coords": coords,
                        "motif": motif,
                        "fwd_seq": fwd_seq,
                        "fwd_tm": fwd_tm,
                        "rev_seq": rev_seq,
                        "rev_tm": rev_tm,
                        "amp_size": amp_size,
                        "score": total_repeat_units
                    })

    # Sort candidates by repeat unit score (highest expansion potential first)
    candidates.sort(key=lambda x: x["score"], reverse=True)

    # Bin into non-overlapping multiplex size ranges
    bins = {
        "150-230bp": [],
        "231-310bp": [],
        "311-450bp": []
    }

    for c in candidates:
        size = c["amp_size"]
        if 150 <= size <= 230:
            bins["150-230bp"].append(c)
        elif 231 <= size <= 310:
            bins["231-310bp"].append(c)
        elif 311 <= size <= 450:
            bins["311-450bp"].append(c)

    # Select top candidates from each size bracket for balanced multiplex panels
    top_panel = []
    for bin_name, bin_list in bins.items():
        selected = bin_list[:10]  # Top 10 per size bracket
        top_panel.extend(selected)
        print(f"  - Size Pool [{bin_name}]: {len(bin_list)} candidate(s) available, selected top {len(selected)}")

    # Write top candidates file
    with open(TOP_MARKERS_FILE, "w") as out:
        out.write("Locus_ID\tChromosome\tGenomic_Coordinates\tSTR_Motif\tExpected_Amplicon_Size\tRepeat_Units_Score\tForward_Primer_5to3\tReverse_Primer_5to3\tForward_Tm\tReverse_Tm\n")
        for p in top_panel:
            out.write(
                f"{p['locus_id']}\t{p['chrom']}\t{p['coords']}\t{p['motif']}\t{p['amp_size']}\t{p['score']}\t"
                f"{p['fwd_seq']}\t{p['rev_seq']}\t{p['fwd_tm']:.2f}\t{p['rev_tm']:.2f}\n"
            )

    print(f"\n✨ Stage 7 Complete: Selected {len(top_panel)} top Tri/Tetra STR diagnostic markers across 3 multiplex size brackets.")
    print(f"📁 Short-listed candidate panel saved to:\n   {TOP_MARKERS_FILE}")


if __name__ == "__main__":
    main()
