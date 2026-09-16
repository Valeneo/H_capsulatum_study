#!/usr/bin/env python3
"""
Step 17: External In-Silico Validation & Specificity Verification Generator
Tools: Python / NCBI Primer-BLAST & UCSC Pre-Configured Validation Suite

Biology & Logic:
- Generates pre-configured validation links and parameters to verify the designed diagnostic primers:
  1. On-Target Verification: Tests primers against NCBI RefSeq database for Histoplasma capsulatum (taxid: 5037).
  2. Host Non-Specific Cross-Reactivity Check: Tests primers against Homo sapiens (taxid: 9606) to guarantee
     zero false-positive amplification from human host clinical samples (blood, tissue biopsy, BAL fluid).
  3. Related Pathogen Specificity Check: Confirms no cross-amplification with other dimorphic fungi
     (Blastomyces, Coccidioides, Aspergillus, Candida).
"""

import os
import urllib.parse
import pandas as pd

KIT_FILE = "results/final_primer_panel/TRUE_final_diagnostic_kit.tsv"
REPORT_FILE = "results/stage6_biophysical_validation/NCBI_UCSC_primer_validation_guide.md"
os.makedirs("results/stage6_biophysical_validation", exist_ok=True)


def main():
    if not os.path.exists(KIT_FILE):
        print(f"❌ Error: Kit file '{KIT_FILE}' not found.")
        return

    df = pd.read_csv(KIT_FILE, sep="\t")

    with open(REPORT_FILE, "w") as out:
        out.write("# Universal STR Diagnostic Panel: External In-Silico Validation Guide\n")
        out.write("**Target Organism**: *Histoplasma capsulatum* (NCBI Taxonomy ID: `5037`)\n\n")
        out.write("This document provides direct validation links and verification parameters to confirm **on-target specificity** and **zero host cross-reactivity** for your thesis/advisor.\n\n")
        out.write("---\n\n")

        for idx, row in df.iterrows():
            locus = row["Locus_ID"]
            motif = row["Motif"]
            fwd = row["Forward_Primer_5to3"]
            rev = row["Reverse_Primer_5to3"]
            fwd_tm = row["Fwd_Tm"]
            rev_tm = row["Rev_Tm"]
            size = row["Expected_Product_Size"]
            alleles = row["Repeat_Alleles_Per_Strain"]

            # 1. NCBI Primer-BLAST On-Target Link (Histoplasma)
            p_on = {
                "PRIMER_LEFT_INPUT": fwd,
                "PRIMER_RIGHT_INPUT": rev,
                "ORGANISM": "Histoplasma capsulatum (taxid:5037)",
                "TOTAL_PRIMER_SPECIFICITY_MISMATCH": 1,
                "MIN_PRODUCT_SIZE": 100,
                "MAX_PRODUCT_SIZE": 1000
            }
            on_target_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?{urllib.parse.urlencode(p_on)}"

            # 2. NCBI Primer-BLAST Host Specificity Link (Human Host)
            p_host = {
                "PRIMER_LEFT_INPUT": fwd,
                "PRIMER_RIGHT_INPUT": rev,
                "ORGANISM": "Homo sapiens (taxid:9606)",
                "TOTAL_PRIMER_SPECIFICITY_MISMATCH": 2,
                "MIN_PRODUCT_SIZE": 100,
                "MAX_PRODUCT_SIZE": 1500
            }
            host_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?{urllib.parse.urlencode(p_host)}"

            # 3. Fungal Cross-Reactivity Check (All Fungi)
            p_fungi = {
                "PRIMER_LEFT_INPUT": fwd,
                "PRIMER_RIGHT_INPUT": rev,
                "ORGANISM": "Fungi (taxid:4751)",
                "TOTAL_PRIMER_SPECIFICITY_MISMATCH": 1,
                "MIN_PRODUCT_SIZE": 100,
                "MAX_PRODUCT_SIZE": 1000
            }
            fungi_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?{urllib.parse.urlencode(p_fungi)}"

            out.write(f"### Locus {idx+1}: `{locus}`\n")
            out.write(f"* **STR Motif**: `{motif}`\n")
            out.write(f"* **Expected Base Product Size**: **{size} bp**\n")
            out.write(f"* **Multi-Clade Allele Profile**: `{alleles}`\n")
            out.write(f"* **Forward Primer (5'->3')**: `{fwd}` (Tm = {fwd_tm}°C)\n")
            out.write(f"* **Reverse Primer (5'->3')**: `{rev}` (Tm = {rev_tm}°C)\n\n")
            out.write(f"#### Direct Validation Actions:\n")
            out.write(f"1. 🎯 **[Click Here to Run On-Target Validation on NCBI Primer-BLAST]({on_target_url})**  \n")
            out.write(f"   *(Verifies single specific band amplification in Histoplasma capsulatum).*\n")
            out.write(f"2. 🛡️ **[Click Here to Check Human Host Cross-Reactivity (Homo sapiens)]({host_url})**  \n")
            out.write(f"   *(Confirms 0 amplicons in human DNA to guarantee no false positives in clinical biopsy/blood).* \n")
            out.write(f"3. 🔬 **[Click Here to Check Fungal Cross-Reactivity (All Fungi)]({fungi_url})**  \n")
            out.write(f"   *(Confirms species-specific diagnosis against other fungal pathogens).*\n\n")
            out.write("---\n\n")

    print("==========================================================================")
    print(f"🎉 VALIDATION GUIDE GENERATED SUCCESSFULLY!")
    print(f"📁 Document saved at: {REPORT_FILE}")
    print("==========================================================================")


if __name__ == "__main__":
    main()
