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

# Atomic lock: mkdir fails if another batch is already working on this sample.
mkdir "$OUT/.lock" 2>/dev/null || { echo "[$RUN] being processed elsewhere, skipping"; exit 0; }
# Disk guard: keep DISK_RESERVE_GB (default 20) free for the system. One sample needs up to ~6 GB of temporary files.
while (( $(df --output=avail -BG /home | tail -1 | tr -dc 0-9) < ${DISK_RESERVE_GB:-20} + 10 )); do
  echo "[$RUN] waiting for disk space"; sleep 120
done
BAM="$OUT/$RUN.bam"
# Trimmed reads may have been deleted to save disk; rebuild them from the .sra if so.
if [[ ! -s "$IN/${RUN}_R1.clean.fastq.gz" || ! -s "$IN/${RUN}_R2.clean.fastq.gz" ]]; then
  echo "[$RUN] regenerating trimmed reads from the .sra"; REGEN=1 bash "$ROOT/scripts/01_qc_trim.sh" "$RUN" "$THREADS"
fi
# BAMs are big and only needed to make the gVCF; delete on exit (keep the small stats files).
trap 'rm -f "$BAM" "$BAM.bai" "$BAM.csi"; rmdir "$OUT/.lock" 2>/dev/null' EXIT

# 1. Align. -R writes the read group: GATK requires every read to say which sample it came from.
#    Then fixmate -m (adds mate-score tags that markdup needs; requires name-sorted input, which bwa gives),
#    sort by position, and markdup to FLAG (not delete) PCR/optical duplicates.
bwa-mem2 mem -t "$THREADS" -R "@RG\tID:$RUN\tSM:$RUN\tPL:ILLUMINA\tLB:$RUN" "$REF" \
    "$IN/${RUN}_R1.clean.fastq.gz" "$IN/${RUN}_R2.clean.fastq.gz" 2>"$OUT/$RUN.bwa.log" \
  | samtools fixmate -m -u - - \
  | samtools sort -u -@ 2 -m 400M -T "$OUT/.sorttmp" - \
  | samtools markdup -@ 2 - "$BAM"
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

# Disk is scarce and the .sra can always regenerate these reads (step 5 does so just in time), so free them now.
if [[ "${KEEP_TRIMMED:-0}" != 1 && -s "$OUT/$RUN.g.vcf.gz.tbi" ]]; then
  rm -f "$IN/${RUN}_R1.clean.fastq.gz" "$IN/${RUN}_R2.clean.fastq.gz"
fi
echo "[$RUN] done"
