# Histoplasma capsulatum Universal STR Diagnostic Study

This repository contains the complete 16-step bioinformatics pipeline, multi-clade whole-genome assemblies, and validated **Universal Multi-Locus STR Diagnostic Kit** for ***Histoplasma capsulatum***.

## 📄 Main Study Report
* **Comprehensive Report**: **[H_capsulatum_study.md](H_capsulatum_study.md)**
* **Master Diagnostic Kit Table**: **[results/final_primer_panel/TRUE_final_diagnostic_kit.tsv](results/final_primer_panel/TRUE_final_diagnostic_kit.tsv)**
* **In-Silico PCR Verification Report**: **[results/stage6_biophysical_validation/complete_insilico_pcr_report.txt](results/stage6_biophysical_validation/complete_insilico_pcr_report.txt)**
* **Biophysical Dimer Screening**: **[results/stage6_biophysical_validation/primer_dimer_screening_report.tsv](results/stage6_biophysical_validation/primer_dimer_screening_report.tsv)**

---

## 🚀 Quick Start & Pipeline Execution

```bash
# 1. Microsatellite Discovery (MISA + TRF)
python3 step1_ssr_discovery.py

# 2. Extract 300bp Consensus Flanks
python3 step2_extract_consensus_flanks.py

# 3. Cross-Genome BLASTn Homology Search
./run_blast_all.sh
python3 step3_filter_conserved.py

# 4. Conserved Sequence Extraction & Primer3 Design
python3 step4_extract_conserved_sequences.py
python3 step5_run_primer3.py
python3 step6_parse_primer3_results.py

# 5. Motif Selection, PIC Scoring & MAFFT MSA
python3 step7_select_top_markers.py
python3 step8_true_polymorphism.py
python3 step9_run_mafft_validation.py

# 6. Truth Table Assembly & Biophysical Validation
python3 step12_get_all_20_repeat_counts.py
python3 step13_assemble_final_truth_table.py
python3 step14_real_pcr_validation.py
python3 step15_complete_pcr_validation.py
python3 step16_primer_dimer_check.py

# 7. Specificity & Host Cross-Reactivity Verification (CLI)
python3 validate_primers_cli.py
```
