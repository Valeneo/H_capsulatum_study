# H. capsulatum comparative genomics pipeline

183 Illumina paired-end WGS isolates (NCBI SRA) -> variants, assemblies, phylogeny, population structure, comparative analysis.
Every stage has a script in `scripts/` and a lesson in `lessons/` explaining how it works and why each choice was made.

Reference: WU24 (GCA_017310585.1). Environments: `conda env hcap` (bwa-mem2, orthofinder, admixture, clipkit); system tools for the rest.

| # | Objective | Scripts | Tools |
|---|---|---|---|
| 1 | Preprocessing / QC | 01_qc_trim.sh, 01b, 01c | vdb-validate, fasterq-dump, fastp |
| 2 | Alignment + variant calling | 00_prepare_reference.sh, 02_align_call.sh, 03_joint_genotype_filter.sh | bwa-mem2, samtools, GATK4 |
| 3 | Assembly + gene prediction | 04_assemble_busco.sh, 05_predict_genes.sh | SPAdes, BUSCO, Augustus |
| 4 | Phylogenomics, population genetics | 06_popgen.sh, 06b_dapc.R | PLINK2, ADMIXTURE, adegenet (DAPC), IQ-TREE2 |
| 5 | Comparative analysis + metadata integration | | OrthoFinder, R |

Verify any finished stage: `python3 scripts/99_verify.py 1|2|3`.
