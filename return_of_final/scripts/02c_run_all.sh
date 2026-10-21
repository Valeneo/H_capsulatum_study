#!/usr/bin/env bash
# Stage 2 for every sample whose stage-1 output exists. JOBS at a time.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
JOBS="${1:-3}"; THREADS="${2:-6}"
mkdir -p "$ROOT/logs"; : > "$ROOT/logs/stage2_failed.txt"

for d in "$ROOT"/results/01_qc/SRR*/; do
  r=$(basename "$d"); [[ -s "$d/${r}_R2.clean.fastq.gz" ]] && echo "$r"
done | xargs -P "$JOBS" -I{} bash -c \
  'bash "$0/scripts/02_align_call.sh" "$1" "$2" >>"$0/logs/stage2_$1.log" 2>&1 || echo "$1" >> "$0/logs/stage2_failed.txt"' \
  "$ROOT" {} "$THREADS"
echo "failed: $(wc -l < "$ROOT/logs/stage2_failed.txt")"
