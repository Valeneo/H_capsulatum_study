#!/usr/bin/env bash
# Objective 2 (cont.): combine all per-sample gVCFs, call jointly, filter to a clean SNP matrix.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$HOME/miniconda3/envs/hcap/bin:$PATH"
REF="$ROOT/reference/WU24.fa"
OUT="$ROOT/results/03_joint"; mkdir -p "$OUT"
THREADS="${1:-8}"

# 1. Sample map: "name<TAB>path" for every finished gVCF.
for f in "$ROOT"/results/02_variants/SRR*/SRR*.g.vcf.gz; do
  printf "%s\t%s\n" "$(basename "$f" .g.vcf.gz)" "$f"
done > "$OUT/sample_map.tsv"
echo "samples: $(wc -l < "$OUT/sample_map.tsv")"

# 2. Joint genotyping, one chromosome at a time (7 independent jobs; bounded memory).
#    GenomicsDBImport stores all gVCFs in an indexed database so GenotypeGVCFs can read them by region.
while read -r CHR; do
  [[ -s "$OUT/$CHR.vcf.gz.tbi" ]] && continue
  rm -rf "$OUT/db_$CHR"
  gatk --java-options "-Xmx10g" GenomicsDBImport \
    --genomicsdb-workspace-path "$OUT/db_$CHR" --sample-name-map "$OUT/sample_map.tsv" \
    -L "$CHR" --reader-threads "$THREADS" --batch-size 50 --tmp-dir "$OUT/tmp" 2>"$OUT/$CHR.import.log"
  gatk --java-options "-Xmx10g" GenotypeGVCFs \
    -R "$REF" -V "gendb://$OUT/db_$CHR" -O "$OUT/$CHR.vcf.gz" -ploidy 1 \
    --include-non-variant-sites false 2>"$OUT/$CHR.genotype.log"
  rm -rf "$OUT/db_$CHR"
done < "$ROOT/reference/nuclear_chromosomes.list"

# 3. Merge chromosomes into one VCF.
bcftools concat -Oz -o "$OUT/all.raw.vcf.gz" "$OUT"/CP*.vcf.gz && bcftools index -t "$OUT/all.raw.vcf.gz"

# 4. Keep biallelic SNPs only (indels are excluded: they are harder to align and to compare across samples).
gatk SelectVariants -R "$REF" -V "$OUT/all.raw.vcf.gz" --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC -O "$OUT/snps.raw.vcf.gz" 2>"$OUT/select.log"

# 5. Site-level hard filters (GATK's own recommended thresholds; we cannot use VQSR - no truth set for this fungus).
#    QD<2 low confidence per unit depth; FS>60/SOR>3 strand bias; MQ<40 poor mapping; the RankSums flag biased alleles.
gatk VariantFiltration -R "$REF" -V "$OUT/snps.raw.vcf.gz" -O "$OUT/snps.flagged.vcf.gz" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRS" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS" 2>"$OUT/filter.log"

# 6. Genotype-level filter: a call made on <5 reads becomes missing (./.) rather than being trusted.
#    Then keep PASS sites present in >=90% of isolates. (A site seen in few isolates cannot be compared.)
bcftools view -f PASS "$OUT/snps.flagged.vcf.gz" \
  | bcftools +setGT -- -t q -n . -i 'FMT/DP<5' \
  | bcftools view -i 'F_MISSING<0.1' -Oz -o "$OUT/snps.final.vcf.gz"
bcftools index -t "$OUT/snps.final.vcf.gz"
echo "final SNPs: $(bcftools view -H "$OUT/snps.final.vcf.gz" | wc -l)"
