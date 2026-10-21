#!/usr/bin/env bash
# Objective 4 (functional characterisation): presence/absence + divergence of virulence, resistance
# and mating-type genes in one isolate's assembly.
# Usage: bash scripts/08_virulence.sh SRR25758805 [threads]
set -euo pipefail

RUN="$1"; THREADS="${2:-4}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ASM="$ROOT/results/04_assembly/$RUN/$RUN.contigs.fa"
Q="$ROOT/reference/virulence_queries.faa"
OUT="$ROOT/results/08_virulence/$RUN"; mkdir -p "$OUT"
[[ -s "$OUT/hits.tsv" ]] && { echo "[$RUN] already done, skipping"; exit 0; }

# tblastn = protein query vs the assembly translated in all six frames.
# We search the ASSEMBLY, not the predicted genes: if Augustus missed or split a gene, the gene is still found.
# A gene reported "absent" therefore really means the DNA is not there (or is too divergent), not "not predicted".
makeblastdb -in "$ASM" -dbtype nucl -out "$OUT/db" -logfile /dev/null
tblastn -query "$Q" -db "$OUT/db" -evalue 1e-10 -num_threads "$THREADS" \
  -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' > "$OUT/raw.tsv"

# Keep the best hit per query and compute how much of the query it covers.
# Coverage matters as much as identity: a paralog can match a small piece of a protein.
sort -k1,1 -k11,11gr "$OUT/raw.tsv" | awk -F'\t' '!seen[$1]++ {printf "%s\t%s\t%.1f\t%.1f\t%s\n",$1,$2,$3,100*($7-$6+1)/$5,$10}' > "$OUT/hits.tsv"
rm -f "$OUT"/db.n*
echo "[$RUN] $(wc -l < "$OUT/hits.tsv") of $(grep -c '>' "$Q") query genes found"
