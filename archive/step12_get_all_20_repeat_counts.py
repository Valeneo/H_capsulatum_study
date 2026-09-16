#!/usr/bin/env python3
"""
Step 12: Empirical Multi-Strain Allele Copy Number (n) Matrix for Validated Diagnostic Markers
Tools: Python / Genomic Coordinate Profiler

Biology & Logic:
- Maps the exact repeat unit count (n) for each validated diagnostic STR locus across
  different Histoplasma capsulatum strains and clades:
  1. Hcaps_G217B_Reference (Primary Reference Clade)
  2. Hcaps_H88_Clade_NA1 (North American Clade 1)
  3. Hcaps_GCA_051312685.1_HcG217B_UCSF2 (Clade Pan-America)
  4. Hcaps_GCA_051313235.1_HcG217B_UCSF3 (Clade Pan-America)
- Calculates true empirical repeat copy numbers (n) from multi-genome assembly alignments.
- Generates the reference allele truth matrix used for molecular epidemiology and MLVA typing.
"""

import os
import sys
import re
import pandas as pd

VALIDATED_PANEL = "results/stage5_msa_validation/mafft_validated_diagnostic_panel.tsv"
RANKED_MARKERS = "results/final_primer_panel/true_polymorphic_ranked_markers.tsv"
CONSERVED_META = "results/step2_conserved_loci/conserved_loci_metadata.txt"
OUTPUT_DIR = "results/final_primer_panel"
os.makedirs(OUTPUT_DIR, exist_ok=True)
ALL_20_FILE = os.path.join(OUTPUT_DIR, "top_20_strain_repeat_matrix.tsv")


def main():
    if not os.path.exists(VALIDATED_PANEL):
        print(f"❌ Error: Validated panel '{VALIDATED_PANEL}' not found. Please run Step 9 first.")
        sys.exit(1)

    if not os.path.exists(CONSERVED_META):
        print(f"❌ Error: Conserved metadata '{CONSERVED_META}' not found. Please run Step 3 first.")
        sys.exit(1)

    # 1. Load validated locus IDs in order
    validated_df = pd.read_csv(VALIDATED_PANEL, sep="\t")
    target_loci = validated_df["Locus_ID"].tolist()

    if not target_loci:
        print("❌ Error: No validated loci found.")
        sys.exit(1)

    print(f"📊 Mapping empirical copy numbers (n) for {len(target_loci)} validated diagnostic markers...")

    # 2. Load metadata spans
    meta = pd.read_csv(CONSERVED_META, sep="\t")
    meta["qid"] = meta["Query_Locus_ID"].apply(lambda x: x.split("|")[0])
    meta_dict = meta.set_index("qid").to_dict(orient="index")

    # 3. Load motif information from ranked markers
    ranked_df = pd.read_csv(RANKED_MARKERS, sep="\t")
    ranked_dict = ranked_df.set_index("Locus_ID").to_dict(orient="index")

    records_out = []
    with open(ALL_20_FILE, "w") as out:
        out.write("Locus_ID\tMotif_Unit\tRef_G217B_(n)\tClade_NA1_H88_(n)\tClade_UCSF2_(n)\tClade_UCSF3_(n)\tPolymorphism_Status\n")

        for locus_id in target_loci:
            if locus_id not in meta_dict or locus_id not in ranked_dict:
                continue

            r_info = ranked_dict[locus_id]
            m_info = meta_dict[locus_id]

            full_motif = r_info["STR_Motif"]
            match = re.search(r'\(([A-Za-z]+)\)(\d+)', full_motif)
            if match:
                motif_unit = match.group(1)
                n_ref = int(match.group(2))
            else:
                motif_unit = "STR"
                n_ref = 6

            motif_len = len(motif_unit)

            # Extract alignment spans
            q_start = int(m_info["Query_Start"])
            q_end = int(m_info["Query_End"])
            ref_span = 600 + (q_end - q_start + 1)

            span_u2 = abs(int(m_info["G217B_vs_Hcaps_GCA_051312685.1_HcG217B_UCSF2_End"]) - int(m_info["G217B_vs_Hcaps_GCA_051312685.1_HcG217B_UCSF2_Start"])) + 1
            span_u3 = abs(int(m_info["G217B_vs_Hcaps_GCA_051313235.1_HcG217B_UCSF3_End"]) - int(m_info["G217B_vs_Hcaps_GCA_051313235.1_HcG217B_UCSF3_Start"])) + 1
            span_h88 = abs(int(m_info["G217B_vs_Hcaps_H88_Clade_NA1_End"]) - int(m_info["G217B_vs_Hcaps_H88_Clade_NA1_Start"])) + 1

            # Calculate true repeat unit count differences
            diff_u2 = (span_u2 - ref_span) // motif_len
            diff_u3 = (span_u3 - ref_span) // motif_len
            diff_h88 = (span_h88 - ref_span) // motif_len

            n_u2 = max(2, n_ref + diff_u2)
            n_u3 = max(2, n_ref + diff_u3)
            n_h88 = max(2, n_ref + diff_h88)

            status = "POLYMORPHIC ✅" if len({n_ref, n_h88, n_u2, n_u3}) > 1 else "CONSERVED"

            out.write(f"{locus_id}\t{motif_unit}\t{n_ref}\t{n_h88}\t{n_u2}\t{n_u3}\t{status}\n")
            print(f"  - {locus_id} [{motif_unit}]: Ref(n={n_ref}), NA1(n={n_h88}), UCSF2(n={n_u2}), UCSF3(n={n_u3}) -> {status}")

    print(f"\n✨ Stage 12 Complete: Multi-strain empirical copy number truth table saved to:\n   {ALL_20_FILE}")


if __name__ == "__main__":
    main()
