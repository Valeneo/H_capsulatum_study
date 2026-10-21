#!/usr/bin/env bash
# Objective 2: alignment and variant calling, one sample.
# clean FASTQ -> bwa-mem2 -> sorted, duplicate-marked BAM -> GATK HaplotypeCaller (haploid) -> gVCF
# Usage: bash scripts/02_align_call.sh SRR25758805 [threads]
set -euo pipefail

RUN="$1"; THREADS="${2:-8}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$HOME/miniconda3/envs/hcap/bin:$PATH"
REF="$ROOT/reference/WU24.fa"
IN="$ROOT/results/01_qc/$RUN"
OUT="$ROOT/results/02_variants/$RUN"
mkdir -p "$OUT"

if [[ -s "$OUT/$RUN.g.vcf.gz.tbi" ]]; then echo "[$RUN] already done, skipping"; exit 0; fi

BAM="$OUT/$RUN.bam"
# BAMs are big and only needed to make the gVCF; delete on exit (keep the small stats files).
trap 'rm -f "$BAM" "$BAM.bai" "$BAM.csi"' EXIT

# 1. Align. -R writes the read group: GATK requires every read to say which sample it came from.
#    Then fixmate -m (adds mate-score tags that markdup needs; requires name-sorted input, which bwa gives),
#    sort by position, and markdup to FLAG (not delete) PCR/optical duplicates.
bwa-mem2 mem -t "$THREADS" -R "@RG\tID:$RUN\tSM:$RUN\tPL:ILLUMINA\tLB:$RUN" "$REF" \
    "$IN/${RUN}_R1.clean.fastq.gz" "$IN/${RUN}_R2.clean.fastq.gz" 2>"$OUT/$RUN.bwa.log" \
  | samtools fixmate -m -u - - \
  | samtools sort -u -@ 4 -T "$OUT/.sorttmp" - \
  | samtools markdup -@ 4 - "$BAM"
samtools index -@ 4 "$BAM"

# 2. Alignment QC: how many reads placed, and how deep we really sequenced the nuclear genome.
samtools flagstat -@ 4 "$BAM" > "$OUT/$RUN.flagstat.txt"
samtools coverage "$BAM" > "$OUT/$RUN.coverage.tsv"

# 3. Variant calling. --ploidy 1: Histoplasma is haploid, so one allele per site, not two.
#    -ERC GVCF: emit a record for EVERY position (variant or "confidently same as reference"),
#    so that later joint genotyping can tell "no variant" apart from "no data" in each isolate.
#    -L: nuclear chromosomes only.
gatk --java-options "-Xmx${HC_XMX:-2g}" HaplotypeCaller \
  -R "$REF" -I "$BAM" -O "$OUT/$RUN.g.vcf.gz" \
  -ploidy 1 -ERC GVCF -L "$ROOT/reference/nuclear_chromosomes.list" \
  >"$OUT/$RUN.gatk.log" 2>&1

echo "[$RUN] done"
