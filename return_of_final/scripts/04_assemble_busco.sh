#!/usr/bin/env bash
# Objective 3: de novo assembly + completeness check, one sample.
# clean FASTQ -> SPAdes -> contig filter -> BUSCO.   No reference is used anywhere in this stage.
# Usage: bash scripts/04_assemble_busco.sh SRR25758805 [threads]
set -euo pipefail

RUN="$1"; THREADS="${2:-8}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IN="$ROOT/results/01_qc/$RUN"
OUT="$ROOT/results/04_assembly/$RUN"
mkdir -p "$OUT"
[[ -s "$OUT/$RUN.contigs.fa" && -s "$OUT/busco_summary.txt" ]] && { echo "[$RUN] already done, skipping"; exit 0; }

# 1. SPAdes builds a de Bruijn graph: reads are cut into k-mers, k-mers that overlap by k-1 are joined,
#    and contigs are paths through the graph (lecture 51-52 assembly theory, but with k-mers instead of whole reads).
#    --careful: runs a mismatch/short-indel corrector afterwards. We do NOT use --isolate because it assumes
#    >100x depth; these samples are ~10-30x.
if [[ ! -s "$OUT/spades/contigs.fasta" ]]; then
  spades.py --careful -1 "$IN/${RUN}_R1.clean.fastq.gz" -2 "$IN/${RUN}_R2.clean.fastq.gz" \
    -t "$THREADS" -m 10 -o "$OUT/spades" >"$OUT/spades.log" 2>&1
fi

# 2. Keep contigs >= 500 bp. Shorter ones are mostly graph debris and break gene prediction.
#    The header is renamed to "<RUN>_contig_N" so every contig name is unique across 183 genomes later.
seqkit seq -m 500 "$OUT/spades/contigs.fasta" | seqkit replace -p '^.*$' -r "${RUN}_contig_{nr}" > "$OUT/$RUN.contigs.fa"
seqkit stats -a -T "$OUT/$RUN.contigs.fa" > "$OUT/assembly_stats.tsv"

# 3. BUSCO: look for the ~4,500 genes that nearly all Onygenales fungi carry exactly once.
#    Found complete & single-copy = the assembly really contains that gene. This is a quality gate
#    (we want >= 95%) and we later reuse these genes as phylogenetic markers.
busco -i "$OUT/$RUN.contigs.fa" -m genome -l onygenales_odb10 -c "$THREADS" \
  -o busco --out_path "$OUT" -f >"$OUT/busco.log" 2>&1
cp "$OUT"/busco/short_summary*.txt "$OUT/busco_summary.txt"

# SPAdes leaves tens of GB of intermediate graphs; drop everything but the contigs.
rm -rf "$OUT/spades" "$OUT/busco/run_onygenales_odb10/hmmer_output" "$OUT/busco/run_onygenales_odb10/miniprot_output"
echo "[$RUN] done"
