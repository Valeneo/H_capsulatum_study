#!/usr/bin/env python3
"""
Master Execution Pipeline for Histoplasma capsulatum STR Diagnostic Kit
Runs all pipeline steps in sequence and logs full terminal output to terminal-output.md
"""

import os
import sys
import subprocess
import time
from datetime import datetime

OUTPUT_MD = "/home/valeneo/Dev/Academia/str/terminal-output.md"

STEPS = [
    ("Step 1: Microsatellite Discovery (MISA + TRF)", "python3 step1_ssr_discovery.py"),
    ("Step 2: Consensus 300bp Flank Extraction", "python3 step2_extract_consensus_flanks.py"),
    ("Step 3A: Cross-Genome BLASTn Alignment", "./run_blast_all.sh"),
    ("Step 3B: Conserved Loci Filtering", "python3 step3_filter_conserved.py"),
    ("Step 4: Conserved Sequence Extraction", "python3 step4_extract_conserved_sequences.py"),
    ("Step 5: Automated Primer3 Design", "python3 step5_run_primer3.py"),
    ("Step 6: Parse Primer3 Output & Pair Assembly", "python3 step6_parse_primer3_results.py"),
    ("Step 7: Tri- & Tetra-Nucleotide Motif Selection", "python3 step7_select_top_markers.py"),
    ("Step 8: Polymorphism Scoring (PIC) & Gel Binning", "python3 step8_true_polymorphism.py"),
    ("Step 9: MAFFT Flank Alignment & Primer Conservation Check", "python3 step9_run_mafft_validation.py"),
    ("Step 12: Extract Repeat Counts Across All Strains", "python3 step12_get_all_20_repeat_counts.py"),
    ("Step 13: Assemble Master Diagnostic Kit Truth Table", "python3 step13_assemble_final_truth_table.py"),
    ("Step 14: Real In-Silico PCR Validation", "python3 step14_real_pcr_validation.py"),
    ("Step 15: Multi-Genome Comprehensive In-Silico PCR Verification", "python3 step15_complete_pcr_validation.py"),
    ("Step 16: Biophysical Primer-Dimer & Hairpin Thermodynamic QC", "python3 step16_primer_dimer_check.py"),
    ("Step 17: External Validation Guide Generator", "python3 step17_external_validation_generator.py"),
    ("CLI Verification: Specificity & Human Host Cross-Reactivity", "python3 validate_primers_cli.py"),
]

def main():
    start_all = time.time()
    header = f"""# Histoplasma capsulatum Universal STR Pipeline - Terminal Output Log

* **Execution Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
* **Target Organism**: *Histoplasma capsulatum* (`txid:5037`)
* **Total Automated Steps**: {len(STEPS)}

---

"""
    # Write initial header
    with open(OUTPUT_MD, "w") as f_out:
        f_out.write(header)

    print("==========================================================================")
    print("🚀 Starting Histoplasma capsulatum STR Pipeline Execution")
    print(f"📝 Logging full terminal output to: {OUTPUT_MD}")
    print("==========================================================================\n")

    for i, (title, cmd) in enumerate(STEPS, 1):
        print(f"[{i}/{len(STEPS)}] Running: {title} ...")
        print(f"    Command: {cmd}")
        step_start = time.time()

        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - step_start

        combined_out = proc.stdout
        if proc.stderr:
            combined_out += ("\n[STDERR]\n" + proc.stderr)

        # Print summary to console
        print(f"    Status: {'✅ SUCCESS' if proc.returncode == 0 else '❌ FAILED'} ({elapsed:.2f}s)\n")

        # Append step output to markdown
        step_section = f"""## {title}

* **Command**: `{cmd}`
* **Status**: `{'SUCCESS (Exit code 0)' if proc.returncode == 0 else f'FAILED (Exit code {proc.returncode})'}`
* **Execution Time**: `{elapsed:.2f} seconds`

```text
{combined_out.strip()}
```

---

"""
        with open(OUTPUT_MD, "a") as f_out:
            f_out.write(step_section)

    total_time = time.time() - start_all
    summary = f"""# 🎉 Pipeline Execution Completed Successfully

* **Total Pipeline Wall Time**: `{total_time:.2f} seconds` (`{total_time/60:.2f} minutes`)
* **Final Master Kit**: `results/final_primer_panel/TRUE_final_diagnostic_kit.tsv`
* **In-Silico PCR Report**: `results/stage6_biophysical_validation/complete_insilico_pcr_report.txt`
* **Primer-Dimer QC Report**: `results/stage6_biophysical_validation/primer_dimer_screening_report.tsv`
* **CLI Specificity Report**: `results/stage6_biophysical_validation/cli_specificity_verification_report.tsv`
"""
    with open(OUTPUT_MD, "a") as f_out:
        f_out.write(summary)

    print("==========================================================================")
    print(f"🎉 Pipeline finished in {total_time:.2f}s ({total_time/60:.2f}m)!")
    print(f"📄 Terminal logs saved to: {OUTPUT_MD}")
    print("==========================================================================")

if __name__ == "__main__":
    main()

