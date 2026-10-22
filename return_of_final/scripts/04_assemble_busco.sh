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

# 0. Disk guard: SPAdes writes 20-40 GB of working files per sample. Wait (do not fail) until there is room.
while (( $(df --output=avail -BG /home | tail -1 | tr -dc 0-9) < ${DISK_RESERVE_GB:-20} + ${SPADES_WORK_GB:-40} )); do
  echo "[$RUN] waiting for disk space (need reserve + SPAdes working space)"; sleep 120
done

# 1. SPAdes builds a de Bruijn graph: reads are cut into k-mers, k-mers that overlap by k-1 are joined,
#    and contigs are paths through the graph (lecture 51-52 assembly theory, but with k-mers instead of whole reads).
#    --careful: runs a mismatch/short-indel corrector afterwards. We do NOT use --isolate because it assumes
#    >100x depth; these samples are ~10-30x.
if [[ ! -s "$OUT/spades/contigs.fasta" ]]; then
  if [[ ! -s "$IN/${RUN}_R1.clean.fastq.gz" || ! -s "$IN/${RUN}_R2.clean.fastq.gz" ]]; then
    echo "[$RUN] regenerating trimmed reads from the .sra (they were deleted to save disk)"
    REGEN=1 bash "$ROOT/scripts/01_qc_trim.sh" "$RUN" "$THREADS"
  fi
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
  --download_path "$ROOT/busco_downloads" --offline -o busco --out_path "$OUT" -f >"$OUT/busco.log" 2>&1
cp "$OUT"/busco/short_summary*.txt "$OUT/busco_summary.txt"

# SPAdes leaves tens of GB of intermediate graphs; drop everything but the contigs.
rm -rf "$OUT/spades" "$OUT/busco/run_onygenales_odb10/hmmer_output" "$OUT/busco/run_onygenales_odb10/miniprot_output"
# The trimmed reads (1-8 GB) are only needed for alignment (step 2) and assembly (this step). If this sample's gVCF exists
# (step 2 done) we no longer need them: delete to keep disk use bounded. The raw .sra can always regenerate them.
# Delete the regenerated trimmed reads only if BOTH the variant calls (step 2) and this assembly verify.
if [[ -s "$ROOT/results/02_variants/$RUN/$RUN.g.vcf.gz.tbi" ]] \
   && python3 "$ROOT/scripts/99_verify.py" 2 --run "$RUN" >/dev/null 2>&1 \
   && python3 "$ROOT/scripts/99_verify.py" 3 --run "$RUN" >/dev/null 2>&1; then
  rm -f "$IN/${RUN}_R1.clean.fastq.gz" "$IN/${RUN}_R2.clean.fastq.gz"; echo "[$RUN] verified; trimmed reads deleted"
else
  echo "$RUN stage3 not verified or step 2 missing; reads kept" >> "$ROOT/logs/needs_review.txt"; echo "[$RUN] reads kept (not verified)"
fi
echo "[$RUN] done"
