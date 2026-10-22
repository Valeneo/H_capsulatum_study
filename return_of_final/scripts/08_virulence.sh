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
tblastn -query "$Q" -db "$OUT/db" -evalue 1e-4 -num_threads "$THREADS" \
  -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' > "$OUT/raw.tsv"

# A gene with introns is split by tblastn into one hit per exon, so the best SINGLE hit covers only part of the protein
# (e.g. a 100%-identical catalase B covered 41% in one hit). Merge all hits from the same contig first, then measure coverage.
python3 - "$OUT/raw.tsv" "$OUT/hits.tsv" <<'PY2'
import sys, collections
raw, out = sys.argv[1], sys.argv[2]
hsp = collections.defaultdict(list); qlen = {}
for l in open(raw):
    q, s, pid, ln, ql, qs, qe, ss, se, ev, bs = l.rstrip("\n").split("\t")
    qlen[q] = int(ql); hsp[(q, s)].append((min(int(qs), int(qe)), max(int(qs), int(qe)), float(pid), int(ln), float(ev)))
best = {}
for (q, s), hs in hsp.items():
    hs.sort(); merged = []
    for a, b, *_ in hs:
        if merged and a <= merged[-1][1] + 1: merged[-1][1] = max(merged[-1][1], b)
        else: merged.append([a, b])
    cov = 100 * sum(b - a + 1 for a, b in merged) / qlen[q]
    ident = sum(h[2] * h[3] for h in hs) / sum(h[3] for h in hs)       # length-weighted identity
    if q not in best or cov > best[q][2]: best[q] = (s, ident, cov, min(h[4] for h in hs))
with open(out, "w") as fh:
    for q, (s, i, c, e) in sorted(best.items()): fh.write(f"{q}\t{s}\t{i:.1f}\t{min(c,100):.1f}\t{e}\n")
PY2
rm -f "$OUT"/db.n*
echo "[$RUN] $(wc -l < "$OUT/hits.tsv") of $(grep -c '>' "$Q") query genes found"
