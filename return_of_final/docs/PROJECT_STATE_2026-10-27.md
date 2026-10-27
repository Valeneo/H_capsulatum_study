# Project record: everything done so far (2026-10-27)

For a session taking over work on the Histoplasma capsulatum comparative-genomics project. Read this together with
`PIPELINE_HANDOFF_BRIEF.md` (what to run) and `../PIPELINE_DESIGN.md` (design and rules). Repo: github.com/Valeneo/H_capsulatum_study, branch main.

## 1. Scientific goal and objectives (v2.0, approved by the supervising professor; the last two were asked to be expanded and have been)
1. Genomic DNA isolation and whole-genome sequencing of H. capsulatum from clinical isolates.
2. Preprocessing for quality control.
3. Alignment and variant calling.
4. Assembly and functional annotation.
5. Phylogenomic analysis, population genetics and functional characterisation: ML phylogeny from single-copy orthologs (IQ-TREE2, model selection, ultrafast bootstrap);
   population structure by PCA (PLINK2), ADMIXTURE and DAPC; virulence genes, antifungal-resistance variants and mating-type idiomorphs (BLAST, HMMER).
6. Comparative analysis integrating genomic evidence with geography and clinical phenotype: gene-family gain/loss between lineages (OrthoFinder orthogroups) and lineage-geography-phenotype concordance.
The analysis uses 183 public Illumina paired-end WGS runs (metadata: `../histoplasma_183_wgs_comprehensive_phenotypic_metadata.csv`; per-isolate clade labels in it are study-level lists, so clades must come from our own tree).
Anchor papers: Teixeira 2016 (global phylogeny) and Sepulveda 2017 (cryptic speciation), PDFs and XML in `../papers/phylogenomics_protocol/`.

## 2. Tool decisions already made (do not reopen without the owner)
- Kept: fastp, Kraken2 (database NOT installed), bwa-mem2, samtools, GATK4, PLINK2, SPAdes, BUSCO (onygenales_odb10, full output), Augustus, OrthoFinder, IQ-TREE2, ADMIXTURE, BLAST, HMMER, adegenet (DAPC).
- Removed on purpose: ASTRAL-III, PhiPack, PopART/pegas, pixy, TreeMix, Dsuite, and the legacy 4-gene MLST harmonisation step. Consequence: single ML tree (no coalescent concordance), no introgression tests.
- Not used: Minimap2 (no long reads in this dataset), Funannotate (needs a GeneMark licence; Augustus with its pre-trained Histoplasma model is called directly), GPU (8 GB card; NVIDIA Parabricks needs 16 GB).
- Wording: the old "dangerous genes" box is called "Virulence Gene Profiling". Flowcharts live in `../plan/` (`histoplasma_simplified.mmd` is the current simple one; `histoplasma_unified.mmd` the old detailed one). Nothing has been pushed to Figma.

## 3. What has been built and verified
Scripts in `../scripts/`, per-step explanations in `../lessons/01-03*.md`. Verifier: `python3 scripts/99_verify.py N [--run SRR]`.
- Step 1 done for all 183 (fastp QC summary: median 59x after trimming, range 13-255x; 9 samples kept 48-69% of reads, still usable).
- Step 2: 7 of 183 gVCFs finished (see the brief). Verified checks: mapping >=90%, depth >=5x, sample name, 7 nuclear chromosomes, haploid genotypes. One sample (SRR13433424, Panama lineage) failed the mapping gate and is flagged, not discarded.
- Step 3 tested on 3 isolates; steps 5, 6, 7 tested on one isolate (SRR25758805): assembly 29.7 Mb, BUSCO 93.1%, 8,308 proteins, 9 of 10 virulence queries found, MAT1-1.
- Corrupted downloads: 7 .sra files failed validation and were replaced and re-validated (originals moved to `../wgs_corrupt/`, kept).

## 4. Findings that change plans
- GATK HaplotypeCaller time follows divergence from the reference, not depth (26x took 211 min, 151x took 483 min, a near-reference 122x took 31 min). Full step 2 is weeks on one slow machine, hence the second-machine handoff and the pending GATK-vs-bcftools question (Li 2023, a C. auris benchmark in the papers folder, found high agreement across 14 pipelines and GATK well balanced).
- Divergent lineages map worse to WU24 (H88 African 80%, G184AR Panama 87%): reference bias, which is why a de novo track exists.
- SPAdes at full depth: roughly an hour per sample including read regeneration and BUSCO, so 6-10 days for 183.
- Disk is scarce: trimmed reads (about 260 GB) are temporary and deleted only after verification; guards keep 20 GB free.

## 5. Open decisions owned by the project owner
GATK vs bcftools for the remaining step 2; a second reference for divergent lineages; deleting `wgs_corrupt/` and eventually the raw `.sra` archives; contamination check for SRR20306115 (55.7% mapped).

## 6. Where things are
Pipeline: `return_of_final/` (README.md has commands and timings). Results and raw data are not in git (`results/`, `../wgs/`). A 140 MB package of the finished step 2 results is at
`return_of_final/handoff/step2_finished_2026-10-27.tgz` on the original machine (copy by USB or cloud).
