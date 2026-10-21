#!/usr/bin/env bash
# Objective 1: preprocessing for quality control.
# One sample at a time: SRA -> FASTQ -> fastp (QC + trim) -> clean FASTQ.gz
# Usage: bash scripts/01_qc_trim.sh SRR25758805 [threads]
set -euo pipefail

RUN="$1"
THREADS="${2:-8}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRA="$ROOT/../wgs/$RUN/$RUN.sra"
OUT="$ROOT/results/01_qc/$RUN"
TMP="$ROOT/results/01_qc/.tmp_$RUN"

# Idempotent: if the clean reads already exist, do nothing. Lets us re-run the
# whole batch after a crash without redoing finished samples.
if [[ -s "$OUT/${RUN}_R1.clean.fastq.gz" && -s "$OUT/${RUN}_R2.clean.fastq.gz" ]]; then
  echo "[$RUN] already done, skipping"; exit 0
fi

mkdir -p "$OUT" "$TMP"
# Raw FASTQ is ~4-5x the size of the .sra, and the disk is small.
# Always delete it when we leave, even on failure.
trap 'rm -rf "$TMP"' EXIT

# 1. Integrity check. A truncated .sra would otherwise fail later with a
#    confusing error, or worse, silently give fewer reads.
vdb-validate "$SRA" >"$OUT/$RUN.vdb-validate.log" 2>&1 \
  || { echo "[$RUN] CORRUPT/INCOMPLETE SRA, see $OUT/$RUN.vdb-validate.log" >&2; exit 1; }

# 2. SRA -> FASTQ. --split-files writes _1/_2 for the two mates of the pair.
fasterq-dump "$SRA" --split-files --threads "$THREADS" --temp "$TMP" -O "$TMP" >/dev/null 2>&1

# 3. fastp: adapter + quality trimming and QC report in one pass.
#   --detect_adapter_for_pe   find adapters from the read overlap, no adapter file needed
#   --cut_right -W 4 -M 20    slide a 4-base window from the 3' end, cut when mean Q < 20
#   -q 20 -u 30               a base is "bad" below Q20; drop a read if >30% of bases are bad
#   -l 50                     drop reads shorter than 50 bp after trimming
#   --trim_poly_g             remove NovaSeq/NextSeq artefact of fake G runs (no-signal = G)
#   NO --dedup                duplicates are marked later at the BAM stage (samtools markdup),
#                             where they are identified by mapping position, which is more accurate
fastp \
  -i "$TMP/${RUN}_1.fastq" -I "$TMP/${RUN}_2.fastq" \
  -o "$OUT/${RUN}_R1.clean.fastq.gz" -O "$OUT/${RUN}_R2.clean.fastq.gz" \
  --detect_adapter_for_pe \
  --cut_right --cut_window_size 4 --cut_mean_quality 20 \
  --qualified_quality_phred 20 --unqualified_percent_limit 30 \
  --length_required 50 \
  --trim_poly_g \
  --thread "$THREADS" \
  --json "$OUT/$RUN.fastp.json" --html "$OUT/$RUN.fastp.html" \
  --report_title "$RUN" 2>"$OUT/$RUN.fastp.log"

echo "[$RUN] done"
