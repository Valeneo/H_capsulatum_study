# Lesson 2: Alignment and variant calling

Objective: *To perform alignment and variant calling.*
Scripts: `00_prepare_reference.sh`, `02_align_call.sh` (one sample), `02c_run_all.sh`, `03_joint_genotype_filter.sh`.
Worked example: real output from SRR25758805.

## The logic

After stage 1 we have millions of short clean reads per isolate, with no idea where on the genome each came from.
To find what makes one isolate different from another, we need a common frame: a **reference genome**.
Place every read on the reference (alignment), then at every position ask "does this isolate carry the reference base or a
different one?" (variant calling). The result is a table of differences on shared coordinates, which is exactly what PCA,
ADMIXTURE and DAPC need later.

## Choosing the reference: WU24

WU24 (GCA_017310585.1) is a *complete* genome (7 chromosomes + mitochondrion, 32.5 Mb, no gaps) of the NAm1 lineage,
the reference lineage in Sepulveda et al. 2017. Complete beats fragmented because reads from a gap-free genome have a
place to land. The trade-off is **reference bias**: isolates far from WU24 (e.g. the African lineage) map worse, and some of
their genes are invisible. That is why Track B (de novo assembly) exists in parallel.

## Step by step

**1. Index the reference (`00_prepare_reference.sh`), once.**
- `bwa-mem2 index`: builds the BWT + FM-index of lectures 19-20, so any read can be searched in the genome without scanning it.
- `samtools faidx`: lets tools jump to a chromosome position directly.
- `gatk CreateSequenceDictionary`: chromosome names, lengths and checksums; GATK refuses to run without it.
- `nuclear_chromosomes.list`: all chromosomes except the mitochondrion. **Why exclude it:** in our test sample the
  mitochondrion has **472x** depth against **12.4x** for the nuclear genome. Hundreds of copies per cell would break any
  depth-based filter and it is not a single-copy nuclear locus.

**2. Align (`bwa-mem2 mem`).** Each read pair is placed on the genome; the pair's two ends must land close together,
facing each other, about one insert length apart (~400 bp here). That pairing is what resolves repeats.
`-R "@RG\tID:..\tSM:.."` writes the **read group**: GATK needs to know which sample each read belongs to.

**3. Fix, sort, mark duplicates (piped, no big temp files).**
- `samtools fixmate -m`: adds mate-score tags that duplicate detection needs.
- `samtools sort`: orders reads by genome position (needed for everything downstream and for indexing).
- `samtools markdup`: *flags* reads that start and end at identical positions (PCR/optical duplicates, lecture 10).
  We flag, not delete: GATK skips flagged reads itself, and nothing is lost if we later change our mind.

**4. Check the alignment before trusting it (SRR25758805):**

| Check | Result | Meaning |
|---|---|---|
| Reads mapped | 99.87% | Same species as the reference, no large contamination |
| Duplicates | 12,495 reads (0.65%) | Library complexity is fine |
| Nuclear depth | 12.4x, every chromosome 11.7-16x | Even coverage |
| Genome covered | 96.7-98.5% per chromosome | The rest is repeats / divergent regions |
| Mitochondrion | 472x | Confirms the exclusion decision |

A sample with <90% mapped reads is a warning sign: wrong species or contamination. `99_verify.py 2` enforces this.

**5. Call variants (`GATK HaplotypeCaller -ploidy 1 -ERC GVCF`).**
- Instead of reading mismatches one position at a time, HaplotypeCaller finds regions that look different, *re-assembles the
  reads there locally* and works out the real sequence. That is why it handles indels and messy regions better than the
  simple mismatch callers in lecture 28.
- `-ploidy 1`: *Histoplasma* is haploid. One allele per site, never 0/1 heterozygotes. Verified: the output contains only
  single-allele genotypes.
- `-ERC GVCF`: writes a record for **every** position, not only variants, saying "confidently matches the reference
  here" or "not enough data here". Later, when 183 isolates are combined, a missing variant means different things:
  "this isolate has the reference base" vs "this isolate was not sequenced there". A normal VCF cannot tell them apart.

**6. Result for SRR25758805:** 36,627 SNPs and ~10,000 indels against WU24 (about 1.1 SNPs per kb), mean depth 12.2 at
variant sites (agrees with the 12.4x alignment depth, so the two checks are consistent). 11% of variant sites have
depth <5; those are turned into "missing" at the filtering step rather than trusted.
The BAM is deleted after calling (it is ~200 MB-1 GB per sample and not needed again).

## Joint genotyping and filtering (`03_joint_genotype_filter.sh`)

1. **GenomicsDBImport + GenotypeGVCFs** call all 183 isolates *together*. A site with weak evidence in one isolate is
   informed by the other 182, and every isolate gets a genotype at every variant site seen in any isolate.
2. **Keep biallelic SNPs only.** Indels are harder to align and to compare across samples.
3. **Site filters** (GATK's recommended hard-filters; VQSR cannot be used because it needs a truth set a fungus does not have):
   QD<2 (low quality per read depth), FS>60 and SOR>3 (strand bias), MQ<40 (reads map ambiguously), plus two rank-sum tests.
4. **Genotype filter:** a call from <5 reads becomes missing. Then keep sites genotyped in >=90% of isolates, since a
   site seen in few isolates cannot be compared across the panel.

The output, `snps.final.vcf.gz`, is the SNP matrix every population-genetics analysis uses.
