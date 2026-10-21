#!/usr/bin/env bash
# Run stage 1 on every sample, JOBS at a time. Failures are listed in logs/stage1_failed.txt
# (a failed integrity check means the .sra must be re-downloaded: scripts/00_redownload.sh)
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
JOBS="${1:-3}"; THREADS="${2:-6}"
mkdir -p "$ROOT/logs"; : > "$ROOT/logs/stage1_failed.txt"

tail -n +2 "$ROOT/config/samples.tsv" | cut -f1 | \
  xargs -P "$JOBS" -I{} bash -c \
  'bash "$0/scripts/01_qc_trim.sh" "$1" "$2" >>"$0/logs/stage1_$1.log" 2>&1 || echo "$1" >> "$0/logs/stage1_failed.txt"' \
  "$ROOT" {} "$THREADS"

python3 "$ROOT/scripts/01b_qc_summary.py"
echo "failed: $(wc -l < "$ROOT/logs/stage1_failed.txt")"
