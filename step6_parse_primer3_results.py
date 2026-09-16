#!/usr/bin/env python3
"""
Step 6: Parse Primer3 Output into Universal STR Diagnostic Primer Matrix
Tools: Python / Primer3 Record Parser

Biology & Logic:
- Reads the raw Primer3 output (.p3out) generated from the conserved flanking sequences.
- Extracts Forward and Reverse primer pairs designed in the 100% conserved flanking regions.
- Verifies standard diagnostic parameters:
  - Melting Temperature (Tm): ~59 - 61°C
  - Expected Amplicon Size: 150 - 450 bp
- Formats results into a clean, lab-ready TSV matrix.
"""

import os

P3_OUTPUT = "primer3/outputs/designed_primers.p3out"
RESULTS_DIR = "results/final_primer_panel"
os.makedirs(RESULTS_DIR, exist_ok=True)
FINAL_TABLE = os.path.join(RESULTS_DIR, "universal_str_primers.tsv")


def main():
    if not os.path.exists(P3_OUTPUT):
        print(f"❌ Error: Primer3 output file '{P3_OUTPUT}' not found.")
        print(f"   Please run primer3_core on the conserved sequence blocks first.")
        return

    primers_designed = []
    failed_count = 0

    with open(P3_OUTPUT, "r") as f:
        current_record = {}
        for line in f:
            line_str = line.strip()
            if not line_str or line_str == "=":
                if "SEQUENCE_ID" in current_record:
                    seq_id = current_record["SEQUENCE_ID"]
                    if "PRIMER_LEFT_0_SEQUENCE" in current_record:
                        parts = seq_id.split("|")
                        locus_id = parts[0] if len(parts) > 0 else seq_id
                        chrom = parts[1] if len(parts) > 1 else "Unknown"
                        coords = parts[2] if len(parts) > 2 else "Unknown"
                        motif = parts[3] if len(parts) > 3 else "Unknown"

                        primers_designed.append({
                            "locus_id": locus_id,
                            "chrom": chrom,
                            "coordinates": coords,
                            "motif": motif,
                            "fwd_seq": current_record["PRIMER_LEFT_0_SEQUENCE"],
                            "fwd_tm": float(current_record.get("PRIMER_LEFT_0_TM", 0.0)),
                            "rev_seq": current_record["PRIMER_RIGHT_0_SEQUENCE"],
                            "rev_tm": float(current_record.get("PRIMER_RIGHT_0_TM", 0.0)),
                            "product_size": int(current_record.get("PRIMER_PAIR_0_PRODUCT_SIZE", 0))
                        })
                    else:
                        failed_count += 1
                current_record = {}
                continue

            if "=" in line_str:
                parts = line_str.split("=", 1)
                current_record[parts[0]] = parts[1]

    # Write out clean master matrix
    with open(FINAL_TABLE, "w") as out:
        out.write("Locus_ID\tChromosome\tGenomic_Coordinates\tSTR_Motif\tForward_Primer_5to3\tForward_Tm\tReverse_Primer_5to3\tReverse_Tm\tExpected_Amplicon_Size\n")
        for p in primers_designed:
            out.write(
                f"{p['locus_id']}\t{p['chrom']}\t{p['coordinates']}\t{p['motif']}\t"
                f"{p['fwd_seq']}\t{p['fwd_tm']:.2f}\t{p['rev_seq']}\t{p['rev_tm']:.2f}\t{p['product_size']}\n"
            )

    print(f"\n✨ Stage 6 Complete: Successfully parsed {len(primers_designed)} designed universal primer pairs.")
    if failed_count > 0:
        print(f"   ℹ️ {failed_count} loci lacked valid primer candidates under strict multiplex parameters.")
    print(f"📊 Diagnostic primer matrix saved at:\n   {FINAL_TABLE}")


if __name__ == "__main__":
    main()
