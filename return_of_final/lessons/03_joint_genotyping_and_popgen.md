# Lesson 3: Joint genotyping, the SNP matrix, and population structure

Scripts: `03_joint_genotype_filter.sh`, `06_popgen.sh`, `06b_dapc.R`.

## Why call all isolates together

Each isolate's gVCF says, at every position, "variant here" or "same as reference here, with confidence N".
Calling isolates one at a time gives 183 separate variant lists that cannot be compared: a site missing from isolate B might
mean "B has the reference base" or "B was not sequenced there". Because the gVCFs recorded *every* position, joint genotyping
can tell the two apart, and a weak signal in one isolate is informed by the other 182.

Mechanically: `GenomicsDBImport` loads all gVCFs into an indexed database (one chromosome at a time to keep memory small),
`GenotypeGVCFs` then writes one genotype per isolate at every site where any isolate differs from the reference.

## The filtering funnel (real numbers, test run on 3 isolates)

| Step | Sites | Why |
|---|---|---|
| Raw variants | 79,804 | everything GATK called |
| Biallelic SNPs | 63,102 | indels are harder to align and compare; sites with 3+ alleles break most population-genetics tools |
| Pass site filters | 58,247 | QD<2, FS>60, SOR>3, MQ<40, rank-sum tests flag low-confidence sites (MQ40 flagged 2,475: reads there map ambiguously, typical of repeats) |
| Final | 32,413 | genotype calls on <5 reads set to missing; keep sites genotyped in >=90% of isolates |

The last step is deliberately strict. With only 3 test isolates, one low-depth isolate removes a site; with 183 isolates the
rule tolerates up to 18 missing, so far more sites survive. Genotypes in the output are only `0` (reference) or `1`
(alternative), never `0/1`, confirming the calling was haploid. 68% of genotypes were alternative: these isolates are
mostly different from the WU24 reference, which is expected for isolates from other lineages.

## From SNP matrix to population structure (`06_popgen.sh`)

1. **Rename chromosomes to 1-7**: ADMIXTURE only accepts numeric chromosome codes. Haploid calls are written as homozygous
   diploid for PLINK (PLINK assumes diploid); this does not invent heterozygosity.
2. **LD pruning** (`--indep-pairwise 50 5 0.2`): drop one of any pair of SNPs with r^2 > 0.2 in a 50-SNP window. Linked SNPs
   count the same evolutionary event many times; left in, one low-recombination region can dominate PC1 and ADMIXTURE
   (which assumes unlinked markers) becomes over-confident.
3. **PCA** (PLINK2): the data's own axes of variation, no assumption about how many groups exist. Look at it first.
4. **ADMIXTURE, K = 1 to 10 with 10-fold cross-validation**: models each isolate as a mix of K ancestral populations.
   The best K has the *lowest* cross-validation error. K is chosen by the data, not assumed.
5. **DAPC** (R, adegenet): `find.clusters` groups isolates without labels; `dapc` then finds the axes that best *separate*
   those groups (PCA finds the axes of most variance, which are not always the ones that separate groups). Cross-validation
   (`xvalDapc`) picks how many principal components to keep so the model does not overfit.

**What to say if asked about clonality:** all three methods assume sufficient recombination to break up linkage. A largely
clonal population has long-range linkage everywhere, so LD pruning removes more than expected and the structure can look
cleaner than it is. The mating-type ratio (MAT1-1 vs MAT1-2, step 7) is the independent evidence for how much recombination
there is.
