# Pipeline design reference (project "str": H. capsulatum comparative genomics)

Location: `/home/valeneo/Dev/Academia/str/return_of_final/`  |  Repo: github.com/Valeneo/H_capsulatum_study (branch main)
Per-step how/why lessons are in `lessons/`; commands and timings in `README.md`. This file is the map and the rules.

## Goal
Phylogenomics and comparative study of 183 Histoplasma capsulatum isolates (Illumina paired-end WGS, public SRA, `../wgs/<SRR>/<SRR>.sra`):
population structure (PCA, ADMIXTURE, DAPC), phylogeny, virulence/mating-type profiling, gene-content comparison, then link to geography and clinical phenotype
(`../histoplasma_183_wgs_comprehensive_phenotypic_metadata.csv`). Reference genome: WU24 (GCA_017310585.1, NAm1, 7 nuclear chromosomes + mitochondrion, 32.5 Mb).

## Two parallel tracks from the same clean reads (the core design idea)
- Track A, reference-based: align to WU24, call haploid variants, joint-genotype -> shared-coordinate SNP matrix (population genetics). Weakness: reference bias (divergent lineages map worse).
- Track B, reference-free: SPAdes assembly, BUSCO, Augustus genes -> gene content and phylogeny. Weakness: no shared coordinates.
Each track covers the other's blind spot.

## Steps (all resumable per sample; scripts in `scripts/`)
| # | What | Script | Done-marker (checked before skipping) |
|---|---|---|---|
| 1 | validate .sra, fasterq-dump, fastp trim + QC | 01_qc_trim.sh, 01c_run_all.sh | results/01_qc/<R>/<R>.fastp.json (written LAST) |
| 2 | bwa-mem2 -> samtools fixmate/sort/markdup -> GATK HaplotypeCaller (-ploidy 1, -ERC GVCF, nuclear chromosomes only) | 02_align_call.sh, 02c_run_all.sh | results/02_variants/<R>/<R>.g.vcf.gz.tbi |
| 3 | joint genotyping, biallelic SNPs, hard filters, DP<5 -> missing, keep sites present in >=90% | 03_joint_genotype_filter.sh | results/03_joint/snps.final.vcf.gz |
| 4 | LD pruning, PCA (PLINK2), ADMIXTURE K=1-10 CV, DAPC (adegenet) | 06_popgen.sh, 06b_dapc.R | results/06_popgen/ |
| 5 | SPAdes (one at a time) + BUSCO onygenales_odb10, full output kept | 04_assemble_busco.sh | results/04_assembly/<R>/busco_summary.txt |
| 6 | Augustus (pre-trained histoplasma_capsulatum model) | 05_predict_genes.sh | results/05_genes/<R>/<R>.faa |
| 7 | tblastn virulence + mating-type queries, merged-exon coverage | 08_virulence.sh, 08b_virulence_summary.py | results/08_virulence/ |
| 8 | BUSCO single-copy genes -> MAFFT -> ClipKIT -> partitioned IQ-TREE2 (outgroup Blastomyces) | 07_*.py/sh | results/07_phylo/ |
`run_pipeline.sh` chains everything; `99_verify.py N [--run SRR]` verifies a stage; `stop.sh` kills safely and removes half-written samples.

## Rules every script follows (keep them when adding code)
1. Idempotent per sample, finished only when the LAST-written file exists. Never test an earlier file.
2. Atomic per-sample lock (`mkdir results/02_variants/<R>/.lock`) so several batches can run safely.
3. Verify before deleting: trimmed reads are deleted only after `99_verify.py` passes; failures stay and are listed in `logs/needs_review.txt`.
4. Disk is scarce, time is not: trimmed reads are temporary (regenerated from the .sra with `REGEN=1 bash scripts/01_qc_trim.sh <R>`). Guards wait if free space < 20 GB reserve + the step's working space. Nothing is depth-capped or subset.
5. Edit running scripts only by atomic replace (write temp file, `mv -f`), never in place: bash reads scripts incrementally.
6. Never `pkill -f` a pattern that appears in your own command line (use the `[x]pattern` bracket trick); use PIDs.
7. Cite where numbers come from: queries for virulence genes come from the WU24 proteome plus NCBI ASU06352.1 (MAT1-2-1); do not invent accessions.

## Shared-machine rules (i7-12700H 20 logical CPUs, 15 GB RAM, btrfs)
- CPU cap 50%: pipeline pinned to logical CPUs 10-19 (`bash cap_cpu.sh 50`; re-run after starting a step from a new shell). Do heavy side work with `taskset -c 0-9`.
- RAM is the tight resource (pipeline holds ~11 GB). No GPU use is possible (8 GB card; Parabricks needs 16 GB).
- btrfs: free space updates lazily; run `sync` before trusting `df`.
- Do not edit `results/` or `../wgs/` by hand.

## Known issues and open decisions (state as of 2026-10-22)
- GATK is far slower than first estimated: 3.5-9.6 h per sample, driven by divergence from the reference, not by depth (26x took 211 min, 151x took 483). Full step 2 is days; whether to switch to bcftools (haploid) or keep GATK is awaiting the user. Plan: compare both on a few finished samples first (Li 2023 C. auris benchmark in `../papers/phylogenomics_protocol/` supports either).
- Low mapping to WU24: SRR20306115 (55.7%, China, GC 38.3%, possible contamination), SRR13433402 (H88, African, 80.4%), SRR13433424 (G184AR, Panama, 86.9%). Divergent lineages show reference bias; the 90% gate flags them for review (logs/needs_review.txt), it does not mean the data is wrong.
- Per-isolate clade labels in the metadata are study-level lists, not per isolate; clades must come from our own tree.
- Not built yet: OrthoFinder comparative analysis, metadata integration and statistics.
- No Kraken2 database (50-70 GB); contamination is judged from mapping rate instead.
