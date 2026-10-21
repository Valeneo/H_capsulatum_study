#!/usr/bin/env bash
# Objective 4 (population genetics): LD-prune the SNP matrix, run PCA (PLINK2) and ADMIXTURE.
# DAPC is in scripts/06b_dapc.R and uses the same pruned SNP set.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$HOME/miniconda3/envs/hcap/bin:$PATH"
OUT="$ROOT/results/06_popgen"; mkdir -p "$OUT"; cd "$OUT"
THREADS="${1:-8}"
VCF="$ROOT/results/03_joint/snps.final.vcf.gz"

# 1. Rename chromosomes CP069114.1 -> 1 ... because ADMIXTURE only accepts numeric chromosome codes,
#    and turn haploid genotypes (1) into diploid-homozygous (1/1). PLINK's data model is diploid;
#    this is a representation trick, it does not invent heterozygosity.
paste <(cat "$ROOT/reference/nuclear_chromosomes.list") <(seq 1 7) > chr_map.txt
bcftools annotate --rename-chrs chr_map.txt "$VCF" -Oz -o snps.numeric.haploid.vcf.gz
bcftools index -t snps.numeric.haploid.vcf.gz
bcftools +fixploidy snps.numeric.haploid.vcf.gz -Oz -o snps.numeric.vcf.gz -- -f 2
bcftools index -t snps.numeric.vcf.gz

# 2. Import. --chr-set 7 no-xy no-mt tells PLINK there are 7 plain chromosomes (it assumes human otherwise).
plink2 --vcf snps.numeric.vcf.gz --double-id --chr-set 7 no-xy no-mt \
  --set-all-var-ids '@:#' --make-pgen --out all

# 3. LD pruning. Windows of 50 SNPs, step 5, drop one of any pair with r^2 > 0.2.
#    Linked SNPs are the same information counted repeatedly; left in, one low-recombination
#    region can dominate PC1 and ADMIXTURE (which assumes unlinked markers) is over-confident.
plink2 --pfile all --chr-set 7 no-xy no-mt --indep-pairwise 50 5 0.2 --out pruned
plink2 --pfile all --chr-set 7 no-xy no-mt --extract pruned.prune.in --make-bed --out pruned_bed
echo "SNPs before/after pruning: $(grep -vc '^#' all.pvar) / $(wc -l < pruned.prune.in)"
# Same pruned sites, still haploid, for DAPC in R.
awk -F: '{print $1"\t"$2}' pruned.prune.in > pruned.sites.tsv
bcftools view -T pruned.sites.tsv snps.numeric.haploid.vcf.gz -Oz -o pruned.haploid.vcf.gz

# 4. PCA, 10 components. Descriptive: no assumption about how many populations exist.
plink2 --bfile pruned_bed --chr-set 7 no-xy no-mt --pca 10 --out pca

# 5. ADMIXTURE for K = 1..10 with 10-fold cross-validation. The best K is the one with the LOWEST CV error.
for K in $(seq 1 10); do
  admixture --cv=10 -j"$THREADS" -s "$K$K" pruned_bed.bed "$K" 2>&1 | tee "admixture_K$K.log" | grep -E '^CV'
done | tee admixture_cv.txt
