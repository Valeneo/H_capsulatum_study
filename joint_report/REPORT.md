# Clade-wise SNP calling and annotation in *Histoplasma capsulatum*: results summary (draft)

Location: `msc/joint/report/` (figures, tables, scripts). Source data are untouched in `msc/joint/`.
The pipeline project (QC, per-sample calling, queue handling) is linked from `msc/str`.

## 1. Headline results

- 180 of 183 isolates were analysed in four clade sets: **NAm1 (n = 117), African (29), Latin American (23), G217B (11)**.
- Each set was called against its own reference: the complete WU24 genome for NAm1, and a de novo assembly for each of the other three. SNPs were then annotated with SnpEff.
- SNP sites per set: NAm1 1,960,029; African 1,344,019; Latin American 1,822,533; G217B 520,860 (all SNPs, no indels in the filtered sets).
- Predicted effects are dominated by non-coding (MODIFIER) classes (91-93%). HIGH-impact effects are 0.05-0.10% of effects: NAm1 8,365; African 4,646; G217B 1,445.
- Missense:silent ratio is below 1 in all annotated sets (NAm1 0.78, African 0.83, G217B 0.66). Ts/Tv is 2.35-2.48.
- **Latin American annotation failed** (no contig matched the SnpEff database), so no effect results exist for those 23 isolates yet.
- **Several sets mix lineages**: within NAm1, 94 isolates carry about 29% non-reference calls while 23 carry about 1%; within G217B, 4 carry about 53% while 7 are identical to the reference.

## 2. What was done (to confirm: items marked [CONFIRM])

1. De novo assembly per clade (African, Latin American, G217B); the NAm1 set uses the complete WU24 genome. [CONFIRM: which isolate or reads were assembled for each clade and the assembler version; the contig names (NODE_x_length_y_cov_z) indicate SPAdes]
2. Reads of each isolate aligned to its clade reference and haploid variants called with GATK, then joint genotyping per clade, filtering to SNPs (`snps.final.vcf.gz`). [CONFIRM: aligner, GATK version, filters and the rule that put each isolate in a set]
3. SnpEff 5.2e annotation with one database per clade. Gene models were transferred with liftoff: the H88 annotation onto the African assembly, and the published G217B annotation onto the G217B assembly; NAm1 uses the WU24 annotation directly. **Latin American has no database:** no annotated genome exists for it, so the run fell back to the WU24 database, whose coordinates do not match the Latin American contig names. (Source: the pipeline session's account; the liftoff outputs are not in this folder, so I could not verify them here.) [CONFIRM: liftoff version and settings, and the fraction of genes lifted per assembly]
4. Summaries here come from the SnpEff `*_snpeff_summary.csv` and `.genes.txt` files and from `bcftools stats` on each annotated VCF, built by `scripts/00_bcftools_stats.sh` and `scripts/01_make_report.R`.

## 3. Results by figure

**Figure 1 (overview).** A: isolates per set. B: divergence of each isolate from its clade reference (non-reference SNPs per kb of reference; one dot per isolate). C: share of predicted effects by impact class, log scale. D: missense and silent effect counts. Latin American is absent from C and D because its annotation failed.

**Figure 2 (genomic context).** A: where SNP effects fall (upstream, downstream, intergenic, UTR, exon, intron, splice). B: coding-region effect classes per 10,000 SNP sites. About 77-79% of effects are upstream or downstream gene variants, which is typical when genes are densely packed. NAm1 has a larger intergenic share (11.9% against 3.4-3.5%) because its annotation has no UTR class; do not compare these two directly.

**Figure 3 (genes).** Genes with the most HIGH-impact variant sites in NAm1, African and G217B. Of the genes listed in each SnpEff gene table, 3,522 of 8,833 (NAm1), 2,658 of 12,025 (African) and 954 of 6,800 (G217B) carry at least one HIGH-impact site. These are candidates only; see section 4.

**Figure 4 (per isolate).** A: divergence from the clade reference. B: missing genotypes (isolates above 15% are labelled). C: mean depth at SNP sites. Median depth is 32x (NAm1), 55x (African), 36x (Latin American), 17x (G217B). Eleven isolates have more than 10% missing genotypes. Six have more than 30%, and five of those six have mean depth of only 5-7x, consistent with low-coverage genotypes being set to missing.

**Figure 5 (reliability).** A: contigs in each reference (8; 16,664; 25,810; 270,166). B: share of SNPs with SnpEff gene-model warnings (NAm1 26.5%, African 27.8%, G217B 8.7%) and the Latin American failure. C: Ts/Tv by set. D: allele frequency spectrum; most sites are rare, and NAm1 has a second peak at 60-80% that matches the group of 94 divergent isolates.

Tables: `tables/clade_summary.csv`, `effects_by_impact.csv`, `effects_by_functional_class.csv`, `effects_by_region.csv`, `effects_by_type.csv`, `top15_high_impact_genes_per_clade.csv`, `gene_level_totals.csv`, `per_isolate_summary.csv`, `divergence_groups_by_clade.csv`, `ts_tv_by_clade.csv`, `isolates_not_in_any_clade_set.csv`.

## 4. Problems to fix or disclose before this is reported

1. **Latin American annotation is invalid.** No annotated genome exists for this lineage, so the annotation fell back to the WU24 database. Its gene coordinates do not match the Latin American assembly's 25,810 contigs, so every variant carries `ERROR_CHROMOSOME_NOT_FOUND`, the missense:silent ratio shows 0, and effects equal variants. The SNP calls themselves are unaffected, and only the effect predictions are void. Options: (a) lift the WU24 or H88 annotation onto the Latin American assembly with liftoff, as was done for African and G217B, build a database and rerun SnpEff on the existing VCF (consistent with the other clades, so my recommendation); (b) predict genes on the assembly with Augustus, which would make this clade's gene models different in kind from the others. A third route is to switch the reference: the NCBI assembly GCA_062029959.1 (*H. suramericanum* MZ5) is contig-level with an N50 near 5 Mb, far less fragmented than the 25,810-contig assembly, but it has no annotation either and would need the calls to be redone.
2. **Three isolates are in no set:** SRR20306115 (China), SRR29905426 (South Africa), SRR29905402 (Egypt). They were among the lowest-mapping isolates in the earlier WU24 and H88 runs (55.7%, 65.6%, 67.2%). The reason for the exclusion is not recorded in these files, so state it explicitly.
3. **Counts are not comparable across clades.** Each set has its own reference of different length (32.5, 42.9, 43.1, 50.2 Mb), contiguity and annotation type, and pooled SNP site counts grow with the number of isolates. Compare per-isolate measures and proportions, not totals.
4. **The de novo references are fragmented and larger than WU24.** G217B has 270,166 contigs and 50.2 Mb; African and Latin American are about 43 Mb. Fragmentation and possible redundant or contaminating contigs can inflate variant counts and break gene models. The assembly sizes should be checked against BUSCO and contig-length statistics.
5. **Gene-model warnings are common.** SnpEff reports incomplete transcripts, missing start or stop codons and internal stop codons for 8.7-27.8% of variants. For African and G217B this is expected, because lifting annotations onto fragmented assemblies truncates genes at contig ends and can introduce frameshifts. Effects on those genes, especially stop-gained and other HIGH-impact calls, are unreliable. Present gene lists as candidates.
6. **Do not cite SnpEff's Hom/Het table.** The genotypes are haploid (0, 1 or missing), so SnpEff labels every non-reference call as "Het". Use non-reference counts per isolate instead.
7. **The clade sets are not pure lineages.** Bimodal divergence (Figure 4A) means NAm1 contains 23 near-identical and 94 divergent isolates, G217B 7 and 4, African 18 and 11, Latin American 2 and 21. A likely reading is that the near-identical group is the lineage of the reference and the divergent group is another lineage, but this is a hypothesis; the phylogeny must confirm it before clade-level conclusions are drawn.
8. **Low-depth isolates have much missing data** (five isolates at 5-7x with 32-60% missing). Treat their calls with caution or exclude them from population-level summaries.

## 5. Suggested wording for the Results section (edit after items 1-4 are resolved)

"Variants were called and annotated separately for four clade sets (NAm1, n = 117; African, n = 29; Latin American, n = 23; G217B, n = 11), each against a clade-specific reference. The sets contained 1.96, 1.34, 1.82 and 0.52 million SNP sites respectively. Across the annotated sets, 91-93% of predicted effects were non-coding, 0.05-0.10% were predicted high impact, and the missense-to-silent ratio was 0.66-0.83. Within-set divergence from the reference was bimodal in all four sets, indicating that each set contains more than one lineage."

## 6. Reproduce

```
bash report/scripts/00_bcftools_stats.sh      # per-sample counts from each annotated VCF (low priority, a few minutes)
Rscript report/scripts/01_make_report.R       # tables and figures
```
