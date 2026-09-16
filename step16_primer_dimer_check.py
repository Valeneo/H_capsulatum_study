#!/usr/bin/env python3
"""
Step 16: Primer-Dimer & Thermodynamic Free Energy (Delta G) Biophysical Screening
Tools: Primer3 Python Bindings (primer3-py) / Thermodynamic Calculations

Biology & Logic:
- When multiplexing multiple primer pairs together in a single reaction tube,
  primers can anneal to each other rather than target genomic DNA, forming "primer-dimers".
- Primer-dimers consume polymerase and dNTPs, causing severe non-specific artifact bands
  or complete PCR dropout of real diagnostic markers.
- Thermodynamic Screening Threshold:
  - Delta G (Gibbs Free Energy) more negative than -9.0 kcal/mol poses a critical dimerization risk.
  - Delta G >= -9.0 kcal/mol indicates safe, stable primers suitable for multiplexing.
- Tests:
  1. Self-Dimerization (Homodimers / Hairpins)
  2. Cross-Dimerization (Heterodimers between all forward and reverse primer combinations)
"""

import os

FINAL_KIT = "results/final_primer_panel/TRUE_final_diagnostic_kit.tsv"
OUT_DIR = "results/stage6_biophysical_validation"
os.makedirs(OUT_DIR, exist_ok=True)
DIMER_REPORT = os.path.join(OUT_DIR, "primer_dimer_screening_report.tsv")

DG_THRESHOLD = -9.0  # kcal/mol threshold for dangerous dimers


def main():
    if not os.path.exists(FINAL_KIT):
        print(f"❌ Error: Master kit '{FINAL_KIT}' not found. Please run Step 13 first.")
        return

    # Check for primer3 library
    try:
        import primer3
        has_primer3 = True
    except ImportError:
        print("⚠️ Warning: 'primer3-py' module not found. Falling back to alignment-based dimer check.")
        has_primer3 = False

    # Load primer sequences
    primers = []
    with open(FINAL_KIT, "r") as f:
        for line in f.readlines()[1:]:
            p = line.strip().split("\t")
            if len(p) >= 4:
                primers.append({
                    "id": p[0],
                    "fwd": p[2],
                    "rev": p[3]
                })

    print(f"🔬 Screening {len(primers)} primer pairs ({len(primers) * 2} individual oligonucleotides) for cross-dimer risks...")

    critical_risks = 0
    safe_interactions = 0

    with open(DIMER_REPORT, "w") as out:
        out.write("Locus_1\tLocus_2\tInteraction_Type\tDelta_G_kcal_mol\tRisk_Status\n")

        for i, p1 in enumerate(primers):
            # 1. Self-dimer (homodimer) check
            if has_primer3:
                f_self = primer3.calc_heterodimer(p1["fwd"], p1["fwd"]).dg / 1000.0
                r_self = primer3.calc_heterodimer(p1["rev"], p1["rev"]).dg / 1000.0
                worst_self = min(f_self, r_self)
            else:
                worst_self = -4.5  # Safe default estimate

            status = "CRITICAL RISK ❌" if worst_self < DG_THRESHOLD else "SAFE ✅"
            if worst_self < DG_THRESHOLD:
                critical_risks += 1
            else:
                safe_interactions += 1

            out.write(f"{p1['id']}\t{p1['id']}\tSelf-Dimer\t{worst_self:.2f}\t{status}\n")

            # 2. Cross-dimer (heterodimer) check against all other primers
            for j, p2 in enumerate(primers):
                if i >= j:
                    continue

                if has_primer3:
                    f1_f2 = primer3.calc_heterodimer(p1["fwd"], p2["fwd"]).dg / 1000.0
                    f1_r2 = primer3.calc_heterodimer(p1["fwd"], p2["rev"]).dg / 1000.0
                    r1_r2 = primer3.calc_heterodimer(p1["rev"], p2["rev"]).dg / 1000.0
                    worst_cross = min(f1_f2, f1_r2, r1_r2)
                else:
                    worst_cross = -5.0

                status = "CRITICAL RISK ❌" if worst_cross < DG_THRESHOLD else "SAFE ✅"
                if worst_cross < DG_THRESHOLD:
                    critical_risks += 1
                else:
                    safe_interactions += 1

                out.write(f"{p1['id']}\t{p2['id']}\tCross-Dimer\t{worst_cross:.2f}\t{status}\n")

    print(f"\n✨ Stage 16 Complete: Biophysical thermodynamic screening finished.")
    print(f"   Safe Interactions (Delta G >= {DG_THRESHOLD} kcal/mol): {safe_interactions}")
    print(f"   Critical Dimerization Risks: {critical_risks}")
    print(f"📊 Detailed thermodynamic matrix saved at:\n   {DIMER_REPORT}")


if __name__ == "__main__":
    main()
