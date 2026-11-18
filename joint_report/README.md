# Clade-wise SNP results (joint calling and SnpEff annotation)

- `REPORT.md`: summary, figure legends, problems to resolve before reporting
- `figures/`: Figures 1-5 (PNG and PDF)
- `tables/`: summary tables behind the figures
- `scripts/`: `00_bcftools_stats.sh` and `01_make_report.R`

The scripts read the annotated VCFs and SnpEff summaries, which are not in this repository (about 1 GB).
They expect to run from the folder that holds those files, with this directory as `report/` inside it.
