#!/usr/bin/env bash
# Objective 3 (cont.): gene prediction on one assembly -> proteome FASTA used by OrthoFinder.
# Usage: bash scripts/05_predict_genes.sh SRR25758805
set -euo pipefail

RUN="$1"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IN="$ROOT/results/04_assembly/$RUN/$RUN.contigs.fa"
OUT="$ROOT/results/05_genes/$RUN"; mkdir -p "$OUT"
[[ -s "$OUT/$RUN.faa" ]] && { echo "[$RUN] already done, skipping"; exit 0; }

# Augustus is a hidden Markov model of what a gene looks like (start codon, exons, introns with GT..AG splice sites,
# stop codon). We use its model PRE-TRAINED on Histoplasma capsulatum, so no RNA-seq training data is needed.
# This is the same predictor Funannotate calls internally; we call it directly because the extra predictors
# Funannotate adds need licences (GeneMark) that are not installed here.
augustus --species=histoplasma_capsulatum --gff3=on --codingseq=on --protein=on "$IN" > "$OUT/$RUN.augustus.gff3"

# Pull the protein sequences out of Augustus's output and give every protein a unique name: <RUN>|<gene>.
perl /usr/share/augustus/scripts/getAnnoFasta.pl --seqfile="$IN" "$OUT/$RUN.augustus.gff3" >/dev/null 2>&1 || true
cat "$OUT"/*.aa | sed "s/^>\(.*\)/>${RUN}|\1/" > "$OUT/$RUN.faa"   # getAnnoFasta chooses its own output name
rm -f "$OUT"/*.aa "$OUT"/*.codingseq "$OUT"/*.cdsexons

echo "[$RUN] genes: $(grep -c '^>' "$OUT/$RUN.faa")"
