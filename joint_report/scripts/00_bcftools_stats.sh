#!/usr/bin/env bash
# Per-sample genotype counts and Ts/Tv from each annotated clade VCF (read-only on the inputs).
set -euo pipefail
export PATH="$HOME/miniconda3/envs/hcap/bin:$PATH"   # falls back to system bcftools if this env is absent
cd "$(dirname "$0")/../.."
for c in G217B LatinAm African NAm1; do
  out="report/tables/raw/${c}_bcftools_stats.txt"; mkdir -p report/tables/raw
  [ -s "$out" ] && continue
  bcftools stats -s - "${c}_joint.ann.vcf.gz" > "$out.tmp" && mv "$out.tmp" "$out"
  echo "$c done"
done
